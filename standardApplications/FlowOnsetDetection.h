/*
 * FlowOnsetDetection.h
 *
 *  Created on: 2013-07-02
 *      Author: ali
 */

#ifndef FLOWONSETDETECTION_H_
#define FLOWONSETDETECTION_H_

#include <cmath>
#include <vector>
#include <stdexcept>

#include "CFL.h"
#include "ConfigXML.h"
#include "StandardAugmentedLagrangian.h"
#include "PrintArguments.h"


template< typename BasicAugmentedLagrangian, typename ForcingRHS >
class FlowOnsetDetection
{
	typedef rheolef::Float Float;
	typedef std::size_t size_t;

public:

	template< typename FieldsPool, typename DirichletBC >
	FlowOnsetDetection( const XMLConfigFile& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		force_rhs(conf.child("rhs_force_term"),fields.Uspace()),
		AL(solver_child(conf),fields,BC),
		algo(solver_child(conf)),
		XML_INIT_VAR(conf,tolerance,"force_param_tolerance"),
		XML_INIT_VAR(conf,initial_dparam,"force_param_delta")
	 {}

	void run()
	{
		Float param = initial_dparam;
		size_t niter = 0;
		while( niter<100 ){
			if( has_flow_for(param) )
				break;
			else
				param += initial_dparam;
			++niter;
		}

		int const n = ceil( log2(initial_dparam/tolerance) );
		Float dparam = initial_dparam/2.;
		param -= dparam;
		for(int i=0; i<n; ++i)
		{
			dparam /= 2.;
			if( has_flow_for(param) )
				param -= dparam;
			else
				param += dparam;
		}
		AL.write_results();
	}

private:

	bool has_flow_for( Float const& param )
	{
		print_args(std::cout,"Param: ",param);
		force_rhs.set_scale_factor(param);
		algo.run(AL,force_rhs);
		rheolef::field const& Gam = AL.get_strainRate_lagrangeMultiplier();
		bool has_flow;
		if( vector_dot(Gam,Gam)==0 ){
			printf(" No flow!\n");
			has_flow = false;
		}
		else {
			printf(" has flow!\n");
			has_flow = true;
		}
		printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n");
		return has_flow;
	}

	XMLConfigFile solver_child( XMLConfigFile const& conf )
	{return conf.child("Solver");}

	ForcingRHS force_rhs;
	BasicAugmentedLagrangian AL;
	standard_augmentedLagrangian_algo algo;

	Float const tolerance;
	Float const initial_dparam;
};


#endif /* FLOWONSETDETECTION_H_ */
