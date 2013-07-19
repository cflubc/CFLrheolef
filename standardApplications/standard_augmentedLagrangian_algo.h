/*
 * standard_augmentedLagrangian_algo.h
 *
 *  Created on: 2013-07-06
 *      Author: ali
 */

#ifndef STANDARD_AUGMENTEDLAGRANGIAN_ALGO_H_
#define STANDARD_AUGMENTEDLAGRANGIAN_ALGO_H_

#include <cstddef>
#include "rheolef/compiler.h"

#include "ConfigXML.h"
#include "ConvergenceMonitor.h"
#include "ResidualTablePrinter.h"

class standard_augmentedLagrangian_algo
{
	typedef rheolef::Float Float;
	typedef std::size_t size_t;

public:

	standard_augmentedLagrangian_algo( const XMLConfigFile& conf ):
		XML_INIT_VAR(conf,convergence_limit,"convergence_limit"),
		residuals_monitor("UTconverge", {"|Un+1-Un|","|Gamdot-Gam|"}),
		XML_INIT_VAR(conf,max_iteration,"max_iteration"),
		n_iterations_without_report( conf("reports_frequency",n_iterations_without_report)-1 ),
		report_header_reprint_frequency( conf.get_if_path_exist({"report_header_reprint_frequency"},30) )
	{}

	template< typename augmentedLag_basic, typename RHSterm >
	void run( augmentedLag_basic& AL, RHSterm const& rhs_manipulator )
	{
		auto program_output = make_residual_table(report_header_reprint_frequency,std::cout,10,16,16);
		size_t niter(0);

		AL.reset_lagrangeMultipliers();
		AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		rhs_manipulator.add_to_rhs( AL.vel_rhs_const_part() );
		do {
			Float Gamres, Ures;
			AL.iterate_ntimes_report_strain_velocity_change(n_iterations_without_report,Gamres,Ures);
			niter += n_iterations_without_report;
			residuals_monitor.add_point(++niter,{Ures,Gamres});

			program_output.print_header_if_needed("\niteration","|Un+1-Un|L2","|Gamdot-Gam|L2");
			program_output.print(niter,Ures,Gamres);
		} while( (niter<max_iteration) && !residuals_monitor.is_converged(convergence_limit) );

		print_solution_convergence_message( niter<max_iteration );
		AL.update_lagrangeMultipliers_clac_strain_rate_multiplier();
	}

	void save_residual_history_to_file() const
	{residuals_monitor.save_to_file();}

	Float const convergence_limit;
	ConvergenceMonitor residuals_monitor;

	size_t const max_iteration;
	size_t const n_iterations_without_report;
	size_t const report_header_reprint_frequency;
};


#endif /* STANDARD_AUGMENTEDLAGRANGIAN_ALGO_H_ */
