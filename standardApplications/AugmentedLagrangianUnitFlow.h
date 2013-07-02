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

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "PrintArguments.h"
#include "CFLSecantMethod.h"
#include "CFLSteadyAnalyser.h"
#include "BorderFluxCalculator.h"
#include "ResidualTablePrinter.h"
#include "AugmentedLagrangian_basic.h"



template< bool isLinear >
class unitflow_iterator;

template<>
class unitflow_iterator<false>
{
public:
	template< typename UnitFlowApp, typename FieldPool >
	unitflow_iterator( UnitFlowApp *const, FieldPool& ) {}

	template< typename UnitFlowApp >
	static void iterate( UnitFlowApp *const app, rheolef::field const& dirichlet_rhs )
	{app->nonlinear_uniflow_iteration(dirichlet_rhs);}
};

template<>
class unitflow_iterator<true>
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

public:
	template< typename UnitFlowApp, typename FieldPool >
	unitflow_iterator( UnitFlowApp *const app, FieldPool& fields ):
		Uh( calc_normalrhs_flow_helper(app,fields) ),
		normalrhs_Uh(fields.Uh()),
		normalrhs_flowrate( app->get_flowrate() )
	{}

	template< typename UnitFlowApp >
	void iterate( UnitFlowApp *const app, field const& dirichlet_rhs )
	{
		auto& predictor = app->predictor;
		auto& AL = app->AL;
		AL.solve_vel_minization();
		Float const flowrate_discripency = predictor.get_target_val()-app->get_flowrate();
		Float const param = flowrate_discripency/normalrhs_flowrate;
		Uh += param*normalrhs_Uh;
		predictor.set_input(param);
		AL.update_lagrangeMultipliers_fast();
		AL.vel_rhs_var_part() = AL.augmented_lagraniang_rhs() + dirichlet_rhs;
	}

private:
	template< typename UnitFlowApp, typename FieldPool >
	static field& calc_normalrhs_flow_helper( UnitFlowApp *const app, FieldPool& fields )
	{
		auto& predictor = app->predictor;
		auto& AL = app->AL;
		app->AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		app->AL.vel_rhs_var_part() = app->rhs_control.get_rhs();
		app->AL.solve_vel_minization();

		return fields.Uh();
	}

	field& Uh;
	field const normalrhs_Uh;
	Float const normalrhs_flowrate;
};




template< typename VelocityMinimizationSolver, typename VelocityRHSManipulator >
class AugmentedLagrangianUnitFlow
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;
	typedef std::size_t size_t;

	enum : bool { useLinearOptimization =  VelocityMinimizationSolver::isLinear &&
                                               VelocityRHSManipulator::isLinear  };
	typedef unitflow_iterator<useLinearOptimization> UnitFlowIterator;
	friend UnitFlowIterator;

public:

	template< typename FieldsPool, typename DirichletBC >
	AugmentedLagrangianUnitFlow( XMLConfigFile const& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		rhs_control(conf.child("unitflow_rhs_controller"),fields.Uspace()),
		flowrate(conf("flowrate_calculation_edge_name"),fields.Uh()),
		AL(conf.child("AugmentedLagrangian"),fields,BC),
		report_header_reprint_frequency( conf.get_if_path_exist({"report_header_reprint_frequency"},30) ),

		LowResolution_conf( conf.child("LowResolution_step") ),
		seq( LowResolution_conf.child("SteadyAnalyser") ),
		predictor( LowResolution_conf.child("Secant") ),
		XML_INIT_VAR(LowResolution_conf,LR_max_iteration,"max_iteration"),
		LR_n_iterations_without_report( XML_VAL(LowResolution_conf,LR_n_iterations_without_report,"reports_frequency")-1 ),

		HighResolution_conf( conf.child("HighResolution_step") ),
		XML_INIT_VAR(HighResolution_conf,HR_max_iteration,"max_iteration"),
		XML_INIT_VAR(HighResolution_conf,HR_converge_limit,"flowrate_convergence_limit"),

		HR_AugLag_conf( HighResolution_conf.child("AugmentedLag") ),
		XML_INIT_VAR(HR_AugLag_conf,HR_AugLag_converge_limit,"convergence_limit"),
		XML_INIT_VAR(HR_AugLag_conf,HR_AugLag_max_iteration,"max_iteration"),
		XML_INIT_VAR(HR_AugLag_conf,HR_AugLag_min_iteration,"min_iteration"),
		HR_AugLag_n_iterations_without_report( HR_AugLag_conf("reports_frequency",HR_AugLag_n_iterations_without_report)-1 ),
		unitflow(this,fields)
	{}


	void run()
	{
		printf("\n-------------------------------------------------------\n"
				 "|                Low Resolution stage                 |\n"
				 "-------------------------------------------------------\n");
		auto output = make_residual_table(report_header_reprint_frequency,std::cout,10,15,17,12);
		ConvergenceMonitor residuals("lowResolution",{"|Un+1-Un|L2","|Gamdot-Gam|L2","ControlParam"});
		size_t niter = 0;

		AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		field const dirichlet_rhs = AL.vel_rhs_const_part();
		while( niter<LR_max_iteration && !seq.sequence_steady_state_reached(residuals[2]) )
		{
			for(size_t i=0; i<LR_n_iterations_without_report; ++i)
				unitflow.iterate(this,dirichlet_rhs);
			niter += LR_n_iterations_without_report;
			Float Ures, Gamres;
			AL.save_strain_velocity();
			unitflow.iterate(this,dirichlet_rhs);
			AL.report_strain_velocity_change(Gamres,Ures);

			residuals.add_point(++niter,{Ures,Gamres,predictor.get_input()});
			output.print_header_if_needed("\niteration","|Un+1-Un|L2","|Gamdot-Gam|L2","Parameter");
			output.print(niter,Ures,Gamres,predictor.get_input());
		}
		residuals.save_to_file();
		printf("\n>>> Unit flow low resolution stage finishde");
		print_solution_convergence_message( niter<LR_max_iteration );


		if( 0<HR_max_iteration ){
			printf("\n-------------------------------------------------------\n"
					 "|                High Resolution stage                |\n"
					 "-------------------------------------------------------\n");
			high_resolution_stage(dirichlet_rhs);
		}

		AL.update_lagrangeMultipliers_clac_strain_rate_multiplier();
		rheolef::odiststream o(AL.geo_name(),"field");
		write_to_diststream(o,"ControlParam", predictor.get_input(),
				              "Flowrate", get_flowrate() );
		AL.write_fields(o);
		o.close();
	}

	field adapt_criteria() const
	{return AL.adapt_criteria();}


private:

	void high_resolution_stage( field const& dirichlet_rhs )
	{
		predictor.set_tolerance_and_Maxiteration(HR_converge_limit,HR_max_iteration);
		predictor.reset();
		ConvergenceMonitor residuals("highRes_augmentedlag",{"|Un+1-Un|L2","|Gamdot-Gam|L2"});

		while( predictor.not_converged_and_have_iterations_left() )
		{
			auto output = make_residual_table(report_header_reprint_frequency,std::cout,10,15,17);
			size_t niter = 0;
			residuals.clear();
			AL.vel_rhs_const_part() = rhs_control.get_rhs( predictor.get_input() ) + dirichlet_rhs;
			AL.reset_lagrangeMultipliers();
			for(;niter<HR_AugLag_max_iteration;)
			{
				Float Gamres, Ures;
				AL.iterate_ntimes_report_strain_velocity_change(HR_AugLag_n_iterations_without_report,Gamres,Ures);
				niter += HR_AugLag_n_iterations_without_report+1;

				residuals.add_point(niter,{Ures,Gamres});
				output.print_header_if_needed("\niteration","|Un+1-Un|L2","|Gamdot-Gam|L2");
				output.print(niter,Ures,Gamres);
				if( HR_AugLag_min_iteration<niter && residuals.is_converged(HR_AugLag_converge_limit) )
					break;
			}

			Float const flux = get_flowrate();
			printf("[Iter %u] Control parameter: %g, Flowrate: %g\n\n",
					predictor.n_iterations_done(), predictor.get_input(), flux);

			predictor.predict_new_input(flux);
		}
		residuals.save_to_file();
	}

	Float get_flowrate() const
	{return flowrate.calc_flux();}

	void nonlinear_uniflow_iteration( field const& dirichlet_rhs )
	{
		predictor.reset();
		// unit flow secant loop
		while( predictor.not_converged_and_have_iterations_left() )
		{
			AL.build_complete_rhs_and_solve_vel_minimization(
							rhs_control.get_rhs( predictor.get_input() )
					 	 	 	 	 	 	 	 	 	 	);
			predictor.predict_new_input( get_flowrate() );
		}
		AL.update_lagrangeMultipliers_fast();
		AL.vel_rhs_const_part() = AL.augmented_lagraniang_rhs() + dirichlet_rhs;
	}


	VelocityRHSManipulator rhs_control;
	BorderFluxCalculator const flowrate;
	AugmentedLagrangian_basic<VelocityMinimizationSolver> AL;
	size_t const report_header_reprint_frequency;

	XMLConfigFile const LowResolution_conf;
	CFLSteadyAnalyser seq;
	CFLSecantMethod predictor;
	size_t const LR_max_iteration;
	size_t const LR_n_iterations_without_report;

	XMLConfigFile const HighResolution_conf;
	size_t const HR_max_iteration;
	Float  const HR_converge_limit;

	XMLConfigFile const HR_AugLag_conf;
	Float  const HR_AugLag_converge_limit;
	size_t const HR_AugLag_max_iteration;
	size_t const HR_AugLag_min_iteration;
	size_t const HR_AugLag_n_iterations_without_report;

	UnitFlowIterator unitflow;
};



#endif /* AUGMENTEDLAGRANGIANUNITFLOW_H_ */
