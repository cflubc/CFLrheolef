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

template< typename VelocityMinimizationSolver, typename VelocityRHSManipulator >
class AugmentedLagrangianUnitFlow
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;
	typedef std::size_t size_t;

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
		HR_AugLag_n_iterations_without_report( HR_AugLag_conf("reports_frequency",HR_AugLag_n_iterations_without_report)-1 )
	{}

	void run()
	{
		printf("\n-------------------------------------------------------\n"
				 "|                Low Resolution stage                 |\n"
				 "-------------------------------------------------------\n");
		auto output = make_residual_table(report_header_reprint_frequency,std::cout,10,15,15,12);
		ConvergenceMonitor residuals("lowResolution",{"|Un+1-Un|L2","|Tn+1-Tn|L2","ControlParam"});
		size_t niter = 0;

		AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		field const dirichlet_rhs = AL.vel_rhs_const_part();
		while( niter<LR_max_iteration && !seq.sequence_steady_state_reached(residuals[2]) )
		{
			for(size_t i=0; i<LR_n_iterations_without_report; ++i)
				uniflow_iteration(dirichlet_rhs);
			niter += LR_n_iterations_without_report;
			Float Ures, Tres;
			AL.save_stress_velocity();
			uniflow_iteration(dirichlet_rhs);
			AL.report_stress_velocity_change(Tres,Ures);

			residuals.add_point(++niter,{Ures,Tres,predictor.get_input()});
			output.print_header_if_needed("\niteration","|Un+1-Un|L2","|Tn+1-Tn|L2","Parameter");
			output.print(niter,Ures,Tres,predictor.get_input());
		}
		residuals.save_to_file();
		printf("\n>>> Unit flow low resolution stage finishde");
		print_solution_convergence_message( niter<LR_max_iteration );



		printf("\n-------------------------------------------------------\n"
				 "|                High Resolution stage                |\n"
				 "-------------------------------------------------------\n");
		predictor.set_tolerance_and_Maxiteration(HR_converge_limit,HR_max_iteration);
		predictor.reset();
		ConvergenceMonitor HR_residuals("highRes_flowrate",{"ControlParam","FlowRate"});
		residuals.rename_and_init("highRes_augmentedlag",{"|Un+1-Un|L2","|Tn+1-Tn|L2"});

		while( predictor.not_converged_and_have_iterations_left() )
		{
			auto output = make_residual_table(report_header_reprint_frequency,std::cout,10,15,15);
			niter = 0;
			residuals.clear();
			AL.vel_rhs_const_part() = rhs_control.get_rhs( predictor.get_input() ) + dirichlet_rhs;
			for(;niter<HR_AugLag_max_iteration;)
			{
				Float Tres, Ures;
				AL.iterate_ntimes_report_stress_velocity_change(HR_AugLag_n_iterations_without_report,Tres,Ures);
				niter += HR_AugLag_n_iterations_without_report+1;

				residuals.add_point(niter,{Ures,Tres});
				output.print_header_if_needed("\niteration","|Un+1-Un|L2","|Tn+1-Tn|L2");
				output.print(niter,Ures,Tres);
				if( HR_AugLag_min_iteration<niter && residuals.is_converged(HR_AugLag_converge_limit) )
					break;
			}

			Float const flux = get_flowrate();
			printf("[Iter %lu] Control parameter: %g, Flowrate: %g\n\n",
					predictor.n_iterations_done(), predictor.get_input(), flux);
			HR_residuals.add_point(predictor.n_iterations_done(),{predictor.get_input(),flux});

			predictor.predict_new_input(flux);
		}
		AL.update_lagrangeMultipliers_clac_strain_rate_multiplier();

		residuals.save_to_file();
		HR_residuals.save_to_file();
		rheolef::odiststream o(AL.geo_name(),"field");
		write_to_diststream(o,"ControlParam", predictor.get_input(),
				              "Flowrate", get_flowrate() );
		AL.write_fields(o);
		o.close();
	}

	field adapt_criteria() const
	{return AL.adapt_criteria();}

private:

	void uniflow_iteration( field const& dirichlet_rhs )
	{
		predictor.reset();
		// unit flow secant loop
		while( predictor.not_converged_and_have_iterations_left() )
		{
			AL.build_complete_rhs_and_solve_vel_minimization(
					rhs_control.get_rhs(predictor.get_input()) );
			predictor.predict_new_input( get_flowrate() );
		}
		AL.update_lagrangeMultipliers_fast();
		AL.vel_rhs_const_part() = AL.augmented_lagraniang_rhs() + dirichlet_rhs;
	}

	Float get_flowrate() const
	{return fabs( flowrate.calc_flux() );}

	VelocityRHSManipulator rhs_control;
	BorderFluxCalculator const flowrate;
	AugmentedLagrangian_basic<VelocityMinimizationSolver> AL;
//	ConvergenceMonitor residuals;
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
};


#endif /* AUGMENTEDLAGRANGIANUNITFLOW_H_ */
