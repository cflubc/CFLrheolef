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

	template< typename Curve >
	T calc_curve_parameter_increment( Curve const& crv, T const& t ) const {
		T const dx = crv.dx(t);
		T const dy = crv.dy(t);
		T const dt = (dx*dx+dy*dy)/fabs(crv.ddy(t)*dx-crv.ddx(t)*dy)*delta_teta_limit;
		return dt;
	}

public:
	ParametricCurveMeshGen( T const& limit = .25 ):
		delta_teta_limit(limit)
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

	std::vector<T> points;
	T t = begin;
	if(include_begin_point)
		points.push_back(begin);

	while( t<end )
	{
		T const dt_predict = calc_curve_parameter_increment(crv,t);
		if( second_order ){
			T const dt_correct = calc_curve_parameter_increment(crv,t+dt_predict);
			t += .5*(dt_predict+dt_correct);
		}
		else
			t += dt_predict;
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

