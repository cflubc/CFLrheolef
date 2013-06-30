/*
 * CFLParametricMeshGen.h
 *
 *  Created on: 2013-06-30
 *      Author: ali
 */

#ifndef CFLPARAMETRICMESHGEN_H_
#define CFLPARAMETRICMESHGEN_H_

#include <cstddef>
#include "rheolef/compiler.h"

#include "ConfigXML.h"
#include "ParametricCurveMeshGen.h"



class CFLConstCurveIntegrator
{
	typedef rheolef::Float Float;
	ConstantCurveIntegrator<Float> I;

public:

	CFLConstCurveIntegrator( XMLConfigFile const& conf ):
		I( conf("npoints_on_curve",std::size_t()))
	{}

	void set_range( Float const& beg, Float const& end )
	{ I.set_range(beg,end); }

	template< typename Curve >
	Float calc_delta_of_curve_parameter( Curve const& crv, Float const& t ) const
	{ return I.calc_delta_of_curve_parameter(crv,t); }

};



class CFLCurveIntegrator
{
	typedef rheolef::Float Float;
	CurveIntegrator<Float,SecondOrderCurveIntegrator> const I;

public:

	CFLCurveIntegrator( XMLConfigFile const& conf ):
		I( conf("max_dTheta",Float()), conf("max_ds",Float()) )
	{}

	void set_range( Float const& beg, Float const& end ) const
	{ I.set_range(beg,end); }

	template< typename Curve >
	Float calc_delta_of_curve_parameter( Curve const& crv, Float const& t ) const
	{ return I.calc_delta_of_curve_parameter(crv,t); }
};


#endif /* CFLPARAMETRICMESHGEN_H_ */
