/*
 * FixedFlowrateIterator.h
 *
 *  Created on: Sep 16, 2013
 *      Author: ali
 */

#ifndef FIXEDFLOWRATEITERATOR_H_
#define FIXEDFLOWRATEITERATOR_H_

#include "rheolef.h"
#include "CFLSecantMethod.h"

template< bool ExploitLinearityOfSystem >
class FixedFlowrateIterator;

template<>
class FixedFlowrateIterator<false>
{
	typedef rheolef::Float Float;
	CFLSecantMethod predictor;
public:
	template< typename UnitFlowApp, typename FieldPool >
	FixedFlowrateIterator( UnitFlowApp *const app, XMLConfigFile const& conf, FieldPool& ):
		predictor( conf.child("Secant"), app->target_flowrate )
	{}

	template< typename UnitFlowApp >
	Float iterate( UnitFlowApp *const app )
	{
		predictor.reset();
		// unit flow secant loop
		while( predictor.not_converged_and_have_iterations_left() )
		{
			app->vel_rhs = app->vel_rhs_const_part + app->rhs_control.get_rhs( predictor.get_input() );
			app->velocity_minimizer.solve(app->vel_rhs);
			predictor.predict_new_input( app->get_flowrate() );
		}
		return predictor.get_input();
	}

	void finalize_iterations( Float const& ) const {}
};


template<>
class FixedFlowrateIterator<true>
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

public:
	template< typename UnitFlowApp, typename FieldPool >
	FixedFlowrateIterator( UnitFlowApp *const app, XMLConfigFile const& conf, FieldPool& fields ):
		Uh( calc_normalrhs_flow_helper(app,fields) ),
		Ph(fields.Ph()),
		normalrhs_Usolution(fields.Uh()),
		normalrhs_Psolution(fields.Ph()),
		normalrhs_flowrate( app->get_flowrate() )
	{}

	template< typename UnitFlowApp >
	Float iterate( UnitFlowApp *const app )
	{
		app->velocity_minimizer.solve( app->vel_rhs_const_part );
		Float const flowrate_discripency = app->target_flowrate - app->get_flowrate();
		Float const param = flowrate_discripency/normalrhs_flowrate;
		Uh += param*normalrhs_Usolution;
		return param;
	}

	void finalize_iterations( Float const& x )
	{Ph += x*normalrhs_Psolution;}

private:
	template< typename UnitFlowApp, typename FieldPool >
	static field& calc_normalrhs_flow_helper( UnitFlowApp *const app, FieldPool& fields )
	{
		app->vel_rhs = app->rhs_control.get_rhs();
		app->velocity_minimizer.solve( app->vel_rhs );
		return fields.Uh();
	}

	field& Uh;
	field& Ph;
	field const normalrhs_Usolution;
	field const normalrhs_Psolution;
	Float const normalrhs_flowrate;
};



#endif /* FIXEDFLOWRATEITERATOR_H_ */
