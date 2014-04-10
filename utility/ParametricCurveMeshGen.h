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
#include <limits>
#include <cstdio>

#include "MathUtility.h"
#include "BiSectionMethod.h"
#include "PrintArguments.h"

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
 *    T estimate_next_point( Curve const&, T const& ) const;
 *
 * Having an integrator scheme, you can pass it to function @c gen_parametric_curve_mesh
 * to generate mesh for curve.
 */


template< typename T >
class ConstantPointCurveIntegrator
{
	std::size_t const Npoints;
	T dt;

public:
	ConstantPointCurveIntegrator( size_t const n ): Npoints(n) {}

	void set_range( T const& beg, T const& end )
	{dt = (end-beg)/(Npoints+1);}

	template< typename Curve >
	T estimate_next_point( Curve const&, T const& t ) const
	{return t+dt;}

	void report() const {
		printf("\n~~~ Constant curve mesh generation report ~~~\n");
		println_args(std::cout," Number of points: ",Npoints);
		printf("~~~~~~\n\n");

	}
};



template< typename T >
struct curve_first_derivative
{
	T const dx;
	T const dy;
	T const dx2plusdy2;
	T const sqrt_dx2plusdy2;

	template< typename Curve >
	curve_first_derivative( Curve const& crv, T const& t ):
		dx( crv.dx(t) ),
		dy( crv.dy(t) ),
		dx2plusdy2( dx*dx+dy*dy ),
		sqrt_dx2plusdy2( std::sqrt(dx2plusdy2) )
	{}

	T angle_of_tangentLine_to_curve( ) const
	{
		T const angle = std::acos( dx/sqrt_dx2plusdy2 );
		return (0<=dy) ? angle : 2.*PI-angle;
	}

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


/** if the curve is straight line, we get a division by zero in
 * @c calc_curve_parameter_increment_2get_const_twist functions.
 * So overload @c estimate_next_point function for the straight line.
 * Here forward declare the class to use it for defining the overload.*/
struct straight_line;

template< typename T, typename IntegrationMethod >
class CurveIntegrator
{
	typedef std::size_t size_t;
public:

	CurveIntegrator( T const& dt, T const& ds ):
		delta_teta_limit(dt),
		points_distance_limit(ds),
		acceptable_tolerance_of_delta_teta(.1*delta_teta_limit),
		n_normal_point(0),
		n_line_point(0),
		n_nonlinear_point(0)
	{}

	void set_range( T const& beg, T const& end ) const {}

	template< typename Curve >
	T estimate_next_point( Curve const& crv, T const& t )
	{
		curve_first_derivative<T> const df(crv,t);
		T const dt = IntegrationMethod::
				calc_curve_parameter_increment_2get_const_twist(crv,t,delta_teta_limit,df);
		T const new_t = t+dt;
		T const curve_tangent_angle = df.angle_of_tangentLine_to_curve();

		if( discrepency_in_delta_teta(crv,delta_teta_limit,curve_tangent_angle,new_t)
			<=acceptable_tolerance_of_delta_teta )
		{
			if( distance_of_points(crv,t,new_t)<=points_distance_limit ){
				++n_normal_point;
				return new_t;
			}
			else {
			// this shows that curve here is almost straight line, use this to find
			// proper dt which approximately gives distance of points as points_distance_limit
				++n_line_point;
				return delta_of_curve_parameter_for_straight_line(df)+t;
			}
		}
		else
		{
		// this is part of curve with sharp twists and the simple linear
		// approximation is not working, use a nonlinear solver
			++n_nonlinear_point;
			return BiSection::solve(t,new_t,7,[&]( T const& x ){
				return discrepency_in_delta_teta(crv,delta_teta_limit,curve_tangent_angle,x);
			});
		}
	}

	T estimate_next_point( straight_line const& crv, T const& t ) const
	{
		static const T dt =
		delta_of_curve_parameter_for_straight_line( curve_first_derivative<T>(crv,t) );
		return t+dt;
	}

	void report() const
	{
		size_t const n_total = n_normal_point+n_line_point+n_nonlinear_point;
		printf("\n~~~ Robust curve mesh generation report ~~~\n");
		println_args(std::cout,"total points generated: ",n_total);
		print_point("normal points: ",n_normal_point,n_total);
		print_point("line points: ",n_line_point,n_total);
		print_point("nonlinear points: ",n_nonlinear_point,n_total);
		printf("~~~~~~\n\n");
	}

private:

	T const delta_teta_limit;
	T const points_distance_limit;
	T const acceptable_tolerance_of_delta_teta;

	size_t n_normal_point;
	size_t n_line_point;
	size_t n_nonlinear_point;

	T delta_of_curve_parameter_for_straight_line( curve_first_derivative<T> const& d ) const
	{return points_distance_limit/d.sqrt_dx2plusdy2;}

	template< typename Curve >
	static T distance_of_points( Curve const& crv, T const& t1, T const& t2 ){
		return sqrt( sqr( crv.x(t2)-crv.x(t1) ) + sqr( crv.y(t2)-crv.y(t1) ) );
	}

	template< typename Curve>
	static T discrepency_in_delta_teta( Curve const& crv, T const& lim,
			                           T const& t1_angle, T const& t2 ){
		curve_first_derivative<T> fd(crv,t2);
		return std::abs(fd.angle_of_tangentLine_to_curve()-t1_angle)-lim;
	}

	template< typename Curve>
	bool point_has_acceptable_delta_teta( Curve const& crv, T const& t1_angle, T const& t2 ) const {
		return discrepency_in_delta_teta(crv,delta_teta_limit,t1_angle,t2)
					   <= acceptable_tolerance_of_delta_teta;
	}

	static T sqr( T const& t )
	{return t*t;};

	static void print_point( char const*const name, double const& n, double const n_total )
	{
		int percent = 100.*n/n_total;
		print_args(std::cout,name,n," (",percent,"%)\n");
	}
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
		t = Ischeme.estimate_next_point(crv,t);
		points.push_back(t);
	} while( t<range.end );
	// in last iteration always one point out of range is inserted
	points.pop_back();
	Ischeme.report();

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


template< typename T, typename Interval >
void gen_straight_line_mesh(
			straight_line const& crv,
			T const& ds,
			std::vector<T> *const X,
			std::vector<T> *const Y,
			Interval)
{
	// for straight line the interval is always 0 to 1
	Interval const range(0,1);
	// here FirstOrderCurveIntegrator and 0 in the constructor are
	// dummy and not important
	CurveIntegrator<T,FirstOrderCurveIntegrator> const I(0,ds);
	gen_parametric_curve_mesh(crv,I,range,X,Y);
}


#endif /* PARAMETRICCURVEMESHGEN_H_ */

