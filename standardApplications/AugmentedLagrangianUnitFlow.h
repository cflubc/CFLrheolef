/*
 * AugmentedLagrangianUnitFlow.h
 *
 *  Created on: 2013-05-17
 *      Author: ali
 */

#ifndef AUGMENTEDLAGRANGIANUNITFLOW_H_
#define AUGMENTEDLAGRANGIANUNITFLOW_H_

#include <cstdlib>
#include "rheolef.h"

#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "CFLSecantMethod.h"
#include "OutputFormatting.h"
#include "BorderFluxCalculator.h"
#include "AugmentedLagrangian_basic.h"




template< typename VelocityMinimizationSolver, typename VelocityRHSManipulator >
class AugmentedLagrangianUnitFlow
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

public:
	template< typename FieldsPool, typename DirichletBC >
	AugmentedLagrangianUnitFlow( XMLConfigFile const& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		rhs_control(conf.child("PBC"),fields.Uh.get_space()),
		flowrate(conf("flowrate_calculation_edge_name"),fields.Uh),
		AL(conf,fields,BC),
		predictor(conf.child("Secant")),
		max_iteration( conf.atoi("max_iteration") ),
		n_iterations_without_report( conf.atoi_if_exist("reports_frequency",10)-1 ),
		Uchange(fields.Uh),
		residuals_monitor("UTCconverge",conf.atof("convergence_limit"),{"|Un+1-Un|","|Tn+1-Tn|","Control"})
	{}

	void run()
	{
		auto program_output = make_column_output(std::cout,6,14,14,11);

		int niter = 0;
		AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		field const dirichlet_rhs = AL.vel_rhs_const_part();

		do {
			for(int i=0; i<n_iterations_without_report; ++i)
				iterate(dirichlet_rhs);
			niter += n_iterations_without_report;

			Uchange.save_field();
			AL.save_stress();
			iterate(dirichlet_rhs);
			Float const Ures = Uchange.calculate_field_change();
			Float const Tres = AL.report_stress_change();
			residuals_monitor.add_point(++niter,{Ures,Tres,predictor.get_input()});
			program_output.print(niter,Ures,Tres,predictor.get_input());
		} while( niter<max_iteration && !residuals_monitor.is_converged() );

		print_solution_convergence_message( residuals_monitor.is_converged() );
		AL.write_results();
		residuals_monitor.save_to_file();
	}

	void iterate( field const& dirichlet_rhs )
	{
		predictor.reset();
		// unit flow secant loop
		while( predictor.not_converged_and_have_iterations_left() )
		{
			AL.build_complete_rhs_and_solve_vel_minimization(
					rhs_control.get_rhs(predictor.get_input()) );
			Float const flux = rheolef::abs( flowrate.calc_flux() );
			predictor.predict_new_input(flux);
		}
		AL.update_lagrangeMultipliers_fast();
		AL.vel_rhs_const_part() = AL.augmented_lagraniang_rhs() + dirichlet_rhs;
	}

private:
	VelocityRHSManipulator rhs_control;
	BorderFluxCalculator const flowrate;
	AugmentedLagrangian_basic<VelocityMinimizationSolver> AL;
	CFLSecantMethod predictor;

	int const max_iteration;
	int const n_iterations_without_report;
	L2norm_calculator Uchange;
	ConvergenceMonitor residuals_monitor;
};


#endif /* AUGMENTEDLAGRANGIANUNITFLOW_H_ */
