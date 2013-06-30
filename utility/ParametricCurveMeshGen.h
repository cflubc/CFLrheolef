/*
 * ParametricCurveMeshGen.h
 *
 *  Created on: 2013-05-24
 *      Author: ali
 */

#ifndef PARAMETRICCURVEMESHGEN_H_
#define PARAMETRICCURVEMESHGEN_H_

#include <cstddef>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

/**
 * Mesh a parametric curve in form {x(t),y(t)}. To mesh curve in range t=[beg,end]
 * we have to use a CurveIntegrator. Currently two types curve integrator available:
 *
 *  - ConstantCurveIntegrator Devides t to N uniformly distributed points in range
 *  - CurveIntegrator More intelligent, tries to choose points such that the change
 *      in curve direction is a constant given value (e.g. .1 radians).
 *      Has first/second order.
 *
 * Required interface for integrator class if you want to create a new one, functions:
 *
 *  - void set_range( T const& beg, T const& end );
 *  - template< typename Curve >
 *    T calc_delta_of_curve_parameter( Curve const&, T const& ) const;
 *
 * Having an integrator scheme, you can pass it to function @c gen_parametric_curve_mesh
 * to generate mesh for curve.
 */


template< typename T >
class ConstantCurveIntegrator
{
	std::size_t const Npoints;
	T dt;

public:
	ConstantCurveIntegrator( size_t const n ): Npoints(n) {}

	void set_range( T const& beg, T const& end )
	{dt = (end-beg)/(Npoints+1);}

	template< typename Curve >
	T calc_delta_of_curve_parameter( Curve const&, T const& ) const
	{return dt;}
};


template< typename T >
struct curve_first_derivative
{
	T const dx;
	T const dy;
	T const dx2plusdy2;

	template< typename Curve >
	curve_first_derivative( Curve const& crv, T const& t ):
		dx( crv.dx(t) ),
		dy( crv.dy(t) ),
		dx2plusdy2( dx*dx+dy*dy )
	{}
};

struct FirstOrderCurveIntegrator
{
	template< typename T, typename Curve >
	static T calc_curve_parameter_increment_2get_const_twist(
			Curve const& crv,
			T const& t,
			T const& delta_teta_limit,
			curve_first_derivative<T> const& fd )
	{
		return fd.dx2plusdy2/fabs(crv.ddy(t)*fd.dx-crv.ddx(t)*fd.dy)*delta_teta_limit;
	}
};


class SecondOrderCurveIntegrator
{
public:
	template< typename T, typename Curve >
	static T calc_curve_parameter_increment_2get_const_twist(
			Curve const& crv,
			T const& t,
			T const& delta_teta_limit,
			curve_first_derivative<T> const& fd )
	{
		T const dt_predict = compute_dt(crv,t,delta_teta_limit);
		T const dt_correct = compute_dt(crv,t+dt_predict,delta_teta_limit);
		return .5*(dt_predict+dt_correct);
	}

private:

	template< typename T, typename Curve >
	static T compute_dt(
			Curve const& crv,
			T const& t,
			T const& delta_teta_limit )
	{
		curve_first_derivative<T> const d(crv,t);
		return FirstOrderCurveIntegrator::
				calc_curve_parameter_increment_2get_const_twist(crv,t,delta_teta_limit,d);
	}
};


template< typename T, typename IntegrationMethod >
class CurveIntegrator
{
public:
	CurveIntegrator( T const& dt, T const& ds ):
		delta_teta_limit(dt),
		points_distance_limit(ds)
	{}

	void set_range( T const& beg, T const& end ) const {}

	template< typename Curve >
	T calc_delta_of_curve_parameter( Curve const& crv, T const& t ) const
	{
		curve_first_derivative<T> const df(crv,t);
		T const dt = IntegrationMethod::
				calc_curve_parameter_increment_2get_const_twist(crv,t,delta_teta_limit,df);

		if( points_distance_limit<distance_of_points(crv,t,t+dt) )
			// this shows that curve here is almost straight line, use this to find
			// proper dt which approximately gives distance of points as points_distance_limit
			return points_distance_limit/sqrt(df.dx2plusdy2);
		else
			return dt;
	}

private:
	T const delta_teta_limit;
	T const points_distance_limit;

	template< typename Curve >
	static T distance_of_points( Curve const& crv, T const& t1, T const& t2 ){
		return sqrt( sqr( crv.x(t2)-crv.x(t1) ) + sqr( crv.y(t2)-crv.y(t1) ) );
	}

	static T sqr( T const& t )
	{return t*t;};
};



template< typename T,
          typename Curve,
          typename Interval,
          typename IntegrationScheme >
void
gen_parametric_curve_mesh(
			Curve const& crv,
			IntegrationScheme& Ischeme,
			Interval const& range,
			std::vector<T> *const X,
			std::vector<T> *const Y ){

	std::size_t const N = X->size();
	if( N!=Y->size() )
		throw std::logic_error("size of X,Y vectors are different");

	std::vector<T> points;
	points.reserve(300);

	if( range.includes_begin() )
		points.push_back(range.beg);

	Ischeme.set_range(range.beg,range.end);
	T t = range.beg;
	do {
		t += Ischeme.calc_delta_of_curve_parameter(crv,t);
		points.push_back(t);
	} while( t<range.end );
	// in last iteration always one point out of range is inserted
	points.pop_back();

	if( range.includes_end() )
		points.push_back(range.end);

	if( range.was_initially_given_reversed )
		std::reverse( points.begin(), points.end() );

	size_t const new_size = points.size()+N;
	X->resize(new_size);
	Y->resize(new_size);

	auto tpoint = points.cbegin();
	for( auto x=X->begin()+N, y=Y->begin()+N;
		 tpoint!=points.cend();
		 ++x, ++y, ++tpoint )
	{
		*x = crv.x( *tpoint );
		*y = crv.y( *tpoint );
	}
}



#endif /* PARAMETRICCURVEMESHGEN_H_ */

