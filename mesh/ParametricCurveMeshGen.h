/*
 * ParametricCurveMeshGen.h
 *
 *  Created on: 2013-05-24
 *      Author: ali
 */

#ifndef PARAMETRICCURVEMESHGEN_H_
#define PARAMETRICCURVEMESHGEN_H_

#include <cmath>
#include <vector>
#include <stdexcept>



struct ParametricMeshOpts
{
	enum: bool { include_ends=true,  exclude_ends=false,
				 include_begin=true, exclude_begin=false,
				 include_end=true,   exclude_end=false,
		         second_order=true,  first_order=false };
};


template< typename T,
          bool include_begin_point,
          bool include_end_point = include_begin_point
         ,bool second_order = true >
class ParametricCurveMeshGen
{
	T const delta_teta_limit;
	T const points_distance_limit;

	template< typename Curve >
	T calc_curve_parameter_increment( Curve const& crv, T const& t, T *const dx2plusdy2 ) const
	{
		T const dx = crv.dx(t);
		T const dy = crv.dy(t);
		*dx2plusdy2 = sqr(dx)+sqr(dy);
		T const dt = *dx2plusdy2/fabs(crv.ddy(t)*dx-crv.ddx(t)*dy)*delta_teta_limit;
		return dt;
	}

	template< typename Curve >
	T distance_of_points( Curve const& crv, T const& t1, T const& t2 ) const {
		return sqrt( sqr( crv.x(t2)-crv.x(t1) ) + sqr( crv.y(t2)-crv.y(t1) ) );
	}

	T sqr( T const& t ) const
	{return t*t;};

public:
	ParametricCurveMeshGen( T const& dt, T const& ds ):
		delta_teta_limit(dt),
		points_distance_limit(ds)
	{}

	template< typename Curve >
	void gen_mesh(
				Curve const& crv,
				T const& begin,
				T const& end,
				std::vector<T> *const X,
				std::vector<T> *const Y ) const;
};


template< typename T,
          bool include_begin_point,
          bool include_end_point,
          bool second_order >
template< typename Curve >
void
ParametricCurveMeshGen<T,include_begin_point,include_end_point,second_order>
::gen_mesh(
			Curve const& crv,
			T const& begin,
			T const& end,
			std::vector<T> *const X,
			std::vector<T> *const Y ) const {

	if( end<begin )
		throw std::logic_error("begin of range should be less than the end of it for parametric mesh generation");

	std::vector<T> points;
	T t = begin;
	if(include_begin_point)
		points.push_back(begin);

	while( t<end )
	{
		T const t_old = t;
		T dx2plusdy;
		T const first_order_dt = calc_curve_parameter_increment(crv,t_old,&dx2plusdy);

		T dt;
		if( second_order ){
			T dx2plusdy2_correct;
			T const dt_correct = calc_curve_parameter_increment(crv,t_old + first_order_dt,&dx2plusdy2_correct);
			dt = .5*(first_order_dt+dt_correct);
		}
		else
			dt = first_order_dt;

		if( points_distance_limit<distance_of_points(crv,t_old,t_old+dt) )
			// this shows that curve here is almost straight line, use this to find
			// proper dt which approximately gives distance of points as max_distance
			dt = points_distance_limit/sqrt(dx2plusdy);

		t += dt;
		points.push_back(t);
	}
	// in last iteration a point with end<=t might be inserted
	if( !(points.back()<end) )
		points.pop_back();

	if(include_end_point)
		points.push_back(end);

	size_t const N = X->size();
	if( N!=Y->size() )
		throw std::logic_error("size of X,Y vectors is different");
	size_t const new_size = points.size()+N;
	X->resize( new_size );
	Y->resize( new_size );

	auto tpoint = points.begin();
	for(size_t i=N; i<new_size; ++i,++tpoint){
		(*X)[i] = crv.x( *tpoint );
		(*Y)[i] = crv.y( *tpoint );
	}
}



#endif /* PARAMETRICCURVEMESHGEN_H_ */

