/*
 * StandardAugmentedLagrangian.h
 *
 *  Created on: 2013-05-02
 *      Author: ali
 */

#ifndef STANDARDAUGMENTEDLAGRANGIAN_H_
#define STANDARDAUGMENTEDLAGRANGIAN_H_

#include "rheolef.h"

#include "ConfigXML.h"
#include "AugmentedLagrangian_basic.h"
#include "standard_augmentedLagrangian_algo.h"


template<
		typename BasicAugmentedLagrangian,
		typename VelocityRHSManipulator >
class StandardAugmentedLagrangian
{
public:

	template< typename FieldsPool, typename DirichletBC >
	StandardAugmentedLagrangian( const XMLConfigFile& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		rhs_manipulator(conf.child("source_term"),fields.Uspace()),
		AL(conf,fields,BC),
		algo(conf)
	{}

	void run() {
		algo.run(AL,rhs_manipulator);
		AL.write_results();
		algo.save_residual_history_to_file();
	}

	rheolef::field adapt_criteria() const
	{return AL.adapt_criteria();}

private:
	VelocityRHSManipulator rhs_manipulator;
	BasicAugmentedLagrangian AL;
	standard_augmentedLagrangian_algo algo;
//	StandardAugmentedLagrangianSolver<BasicAugmentedLagrangian> sAL;
};


#endif /* STANDARDAUGMENTEDLAGRANGIAN_H_ */
