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
#include "StandardAugmentedLagrangianSolver.h"


template< typename VelocityMinimizationSolver, typename VelocityRHSManipulator >
class StandardAugmentedLagrangian
{
	typedef rheolef::field field;
public:

	template< typename FieldsPool, typename DirichletBC >
	StandardAugmentedLagrangian( const XMLConfigFile& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		rhs_manipulator(conf.child("source_term"),fields.Uspace()),
		sAL(conf,fields,BC)
	{}

	void run() {
		sAL.solve(rhs_manipulator);
		sAL.write_results();
	}

	field adapt_criteria() const
	{return sAL.adapt_criteria();}

private:
	VelocityRHSManipulator rhs_manipulator;
	StandardAugmentedLagrangianSolver<VelocityMinimizationSolver> sAL;
};


#endif /* STANDARDAUGMENTEDLAGRANGIAN_H_ */
