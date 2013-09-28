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



class CFLConstCurveIntegrator : public ConstantPointCurveIntegrator<rheolef::Float>
{
public:
	CFLConstCurveIntegrator( XMLConfigFile const& conf ):
		ConstantPointCurveIntegrator<rheolef::Float>(
				conf("npoints_on_curve",std::size_t()) )
	{}
};


class CFLCurveIntegrator : public CurveIntegrator<rheolef::Float,SecondOrderCurveIntegrator>
{
	typedef rheolef::Float Float;
public:

	CFLCurveIntegrator( XMLConfigFile const& conf ):
		CurveIntegrator<Float,SecondOrderCurveIntegrator>(
				conf("max_dTheta",Float()), conf("max_ds",Float()) )
	{}
};


#endif /* CFLPARAMETRICMESHGEN_H_ */
