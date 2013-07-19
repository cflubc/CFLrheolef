/*
 * StandardAugmentedLagrangianSolver.h
 *
 *  Created on: 2013-07-06
 *      Author: ali
 */

#ifndef STANDARDAUGMENTEDLAGRANGIANSOLVER_H_
#define STANDARDAUGMENTEDLAGRANGIANSOLVER_H_

#include "rheolef.h"

#include "AugmentedLagrangian_basic.h"
#include "standard_augmentedLagrangian_algo.h"


template< typename VelocityMinimizationSolver >
class StandardAugmentedLagrangianSolver
{
	typedef rheolef::field field;

public:

	template< typename FieldsPool, typename DirichletBC >
	StandardAugmentedLagrangianSolver( const XMLConfigFile& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		AL(conf.child("AugmentedLagrangian"),fields,BC),
		algo(conf)
	{}

	template< typename RHSterm >
	void solve( RHSterm& rhs )
	{algo.run(AL,rhs);}

	void write_results() {
		AL.write_results();
		algo.save_residual_history_to_file();
	}

	field adapt_criteria() const
	{return AL.adapt_criteria();}

	field const& get_strainRate_lagrangeMultiplier() const
	{return AL.get_strainRate_lagrangeMultiplier();}

	AugmentedLagrangian_basic<VelocityMinimizationSolver> AL;
	standard_augmentedLagrangian_algo algo;
};


#endif /* STANDARDAUGMENTEDLAGRANGIANSOLVER_H_ */
