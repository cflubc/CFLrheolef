/*
 * ParametricCurves.h
 *
 *  Created on: 2013-05-26
 *      Author: ali
 *  to avoid name conflicst, precede the name of classes with "shape_"
 */

#ifndef PARAMETRICCURVES_H_
#define PARAMETRICCURVES_H_


#include <cmath>
#include "MathUtility.h"



#define SHAPE_FUNC(fcn,return_expr) \
	double fcn( double const& t ) const \
	{ return return_expr; }


struct shape_ellipse
{
	double const rx;
	double const ry;
	double const xc;
	double const yc;

	shape_ellipse( double const& _a, double const& _b,
				   double const& _x, double const& _y ):
	   rx(_a),
	   ry(_b),
	   xc(_x),
	   yc(_y)
	{}

	shape_ellipse( double const& _a, double const& _b ):
		shape_ellipse(_a,_b,0.,0.)
	{}

	SHAPE_FUNC( x, xc+rx*std::sin(t) )
	SHAPE_FUNC( y, yc+ry*std::cos(t) )

	SHAPE_FUNC( dx, rx*std::cos(t) )
	SHAPE_FUNC( dy,-ry*std::sin(t) )

	SHAPE_FUNC( ddx,-rx*std::sin(t) )
	SHAPE_FUNC( ddy,-ry*std::cos(t) )
};


struct wavy_wall
{
	double const amplitude;
	double const half_of_wavelength;

	wavy_wall( double const& a, double const& wave_length ):
		amplitude(a),
		half_of_wavelength(.5*wave_length)
	{}

	double scale() const
	{return PI/half_of_wavelength;}

	double scaled_t( double const& t ) const
	{return scale()*t;}

	SHAPE_FUNC( x, t )
	SHAPE_FUNC( y, 1.+amplitude*(1.+std::cos(scaled_t(t))) )

	SHAPE_FUNC( dx, 1. )
	SHAPE_FUNC( dy, -amplitude*scale()*std::sin(scaled_t(t)) )

	SHAPE_FUNC( ddx, 0. )
	SHAPE_FUNC( ddy, -amplitude*scale()*scale()*std::cos(scaled_t(t)) )
};


/** Curve parameter of line varies from [0 1] from point 1 to 2 */
struct straight_line
{
	double const x1;
	double const y1;

	double const x2;
	double const y2;

	double const x21;
	double const y21;

	straight_line( double const& _x1, double const& _y1,
			       double const& _x2, double const& _y2 ):
	   x1(_x1),
	   y1(_y1),
	   x2(_x2),
	   y2(_y2),
	   x21(_x2-_x1),
	   y21(_y2-_y1)
	{}

	SHAPE_FUNC( x, x1+t*x21 )
	SHAPE_FUNC( y, y1+t*y21 )

	SHAPE_FUNC( dx, x21 )
	SHAPE_FUNC( dy, y21 )

	SHAPE_FUNC( ddx, 0. )
	SHAPE_FUNC( ddy, 0. )
};


#undef SHAPE_FUNC
#endif /* PARAMETRICCURVES_H_ */
