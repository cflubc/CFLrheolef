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



#define SHAPE_FUNC(fcn,return_expr) \
	double fcn( double const& t ) const \
	{ return return_expr; }


struct shape_ellipse
{
	double const rx;
	double const ry;
	double const xc;
	double const yc;

public:
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

	SHAPE_FUNC( x,xc+rx*std::sin(t) )
	SHAPE_FUNC( y,yc+ry*std::cos(t) )

	SHAPE_FUNC( dx, rx*std::cos(t) )
	SHAPE_FUNC( dy,-ry*std::sin(t) )

	SHAPE_FUNC( ddx,-rx*std::sin(t) )
	SHAPE_FUNC( ddy,-ry*std::cos(t) )
};


#undef SHAPE_FUNC
#endif /* PARAMETRICCURVES_H_ */
