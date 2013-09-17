/*
 * AugmentedLagrangianUnitFlow.h
 *
 *  Created on: 2013-05-17
 *      Author: ali
 */

#ifndef AUGMENTEDLAGRANGIANUNITFLOW_H_
#define AUGMENTEDLAGRANGIANUNITFLOW_H_

#include <cstddef>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "BorderFluxCalculator.h"
#include "ResidualTablePrinter.h"
#include "AugmentedLagrangian_basic.h"
#include "ConvergenceMonitor.h"
#include "BlockSystem_abtb.h"
#include "IncompressibleStokesSolver.h"
#include "GenericIteration.h"
#include "FixedFlowrateIterator.h"



template< typename BasicAugmentedLagrangian,
		  typename VelocityRHSManipulator >
class AugmentedLagrangianUnitFlow
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;
	typedef std::size_t size_t;

	typedef FixedFlowrateIterator<VelocityRHSManipulator::isLinear> UnitFlowIterator;
	friend UnitFlowIterator;

public:

	template< typename FieldsPool, typename DirichletBC >
	AugmentedLagrangianUnitFlow( XMLConfigFile const& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		XML_INIT_VAR(conf,target_flowrate,"target_flowrate"),
		AL(conf,fields,BC),
		velocity_minimizer(conf,fields,BC,AL.augmentation_coef()),
		rhs_control(conf.child("unitflow_rhs_controller"),fields.Uspace()),
		vel_rhs_const_part(fields.Uspace(), 0.),
		vel_rhs(vel_rhs_const_part),
		dirichlet_rhs( init_dirichlet_rhs(fields.Uspace()) ),
		flowrate(conf("flowrate_calculation_edge_name"),fields.Uh()),
		deltaU(fields.Uh()),
		loop(conf),
		base_name( fields.geo_name() ),
		output(30,std::cout,10,15,17,12),
		residuals("lowResolution",{"|Un+1-Un|L2","|Gamdot-Gam|L2","ControlParam"}),
		unitflow(this,conf,fields)
	{}

	void iterate(){
		flowrate_control_param = unitflow.iterate(this);
		AL.update_lagrangeMultipliers_fast();
		vel_rhs_const_part = AL.augmented_lagraniang_rhs() + dirichlet_rhs;
	}

	void iterate_report( size_t const niter, Float& res ){
		deltaU.save_field();
		AL.save_strain();
		iterate();
		Float const Gamres = AL.strain_change();
		Float const Ures = deltaU.calculate_field_change();
		res = rheolef::max(Ures,Gamres);

		output.print_header_if_needed("\niteration","|Un+1-Un|L2","|Gamdot-Gam|L2","Parameter");
		output.print(niter,Ures,Gamres,flowrate_control_param);
		residuals.add_point(niter,{Ures,Gamres,flowrate_control_param});
	}

	void run()
	{
		printf("\n-------------------------------------------------------\n"
				 "|               Fixed flowrate iteration              |\n"
				 "-------------------------------------------------------\n");
		loop(*this);
		unitflow.finalize_iterations(flowrate_control_param);
		AL.update_lagrangeMultipliers_clac_strain_rate_multiplier();
		residuals.save_to_file();
		rheolef::odiststream o(AL.geo_name(),"field");
		write_to_diststream(o,"ControlParam", flowrate_control_param,
				              "Flowrate", get_flowrate() );
		velocity_minimizer.write_results(o);
		AL.write_fields(o);
		o.close();
	}

	field adapt_criteria() const
	{return AL.adapt_criteria();}

private:

	Float get_flowrate() const
	{return flowrate.calc_flux();}

	field init_dirichlet_rhs( rheolef::space const& X ){
		field f(X,0.);
		velocity_minimizer.set_discrete_dirichlet_rhs(f);
		return f;
	}


	Float flowrate_control_param;
	Float const target_flowrate;
	BasicAugmentedLagrangian AL;
	IncompLinearDiffusionStokesSolver<BlockSystem_abtb> velocity_minimizer;
	VelocityRHSManipulator rhs_control;
	field vel_rhs_const_part;
	field vel_rhs;
	field const dirichlet_rhs;
	BorderFluxCalculator const flowrate;
	L2norm_calculator deltaU;

	GenericIteration loop;
	std::string const base_name;
	ResidualTablePrinter<4,std::ostream> output;
	ConvergenceMonitor residuals;
	UnitFlowIterator unitflow;
};



#endif /* AUGMENTEDLAGRANGIANUNITFLOW_H_ */
