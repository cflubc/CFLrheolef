/*
 * StandardAugmentedLagrangian.h
 *
 *  Created on: 2013-05-02
 *      Author: ali
 */

#ifndef STANDARDAUGMENTEDLAGRANGIAN_H_
#define STANDARDAUGMENTEDLAGRANGIAN_H_

#include <iostream>
#include <cstdlib>
#include "rheolef.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "OutputFormatting.h"
#include "ConvergenceMonitor.h"
#include "NormalStressBC_RHS.h"
#include "AugmentedLagrangian_basic.h"



//struct VoidRHS
//{
//	template< typename FieldsPool >
//	VoidRHS( const XMLConfigFile&, FieldsPool& ) {}
//	void add_to_rhs( rheolef::field& ) const {}
//};
struct VoidRHS
{
	VoidRHS( XMLConfigFile const&, rheolef::space const& ) {}
	void add_to_rhs( rheolef::field& ) const {}
};


template< typename VelocityMinimizationSolver, typename VelocityRHSManipulator >
class StandardAugmentedLagrangian
{
	typedef rheolef::Float Float;
	typedef rheolef::field field;

public:
	template< typename FieldsPool, typename DirichletBC >
	StandardAugmentedLagrangian( const XMLConfigFile& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		rhs_manipulator(conf.child("source_term"),fields.Uh.get_space()),
		AL(conf.child("AugmentedLagrangian"),fields,BC),
		max_iteration( conf.atoi("max_iteration") ),
		n_iterations_without_report( conf.atoi_if_exist("reports_frequency",10)-1 ),
		Uchange(fields.Uh),
		residuals_monitor("UTconverge", conf.atof("convergence_limit"), {"|Un+1-Un|","|Tn+1-Tn|"}),
		time_to_print_header( conf.atoi_if_exist("report_header_reprint_frequency",30) )
	{}

	void run()
	{
		auto program_output( make_column_output(std::cout,10,16,16) );
		int niter(0);

		AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		rhs_manipulator.add_to_rhs( AL.vel_rhs_const_part() );
		do {
			// iterations without reporting
			AL.iterate_ntimes(n_iterations_without_report);
			niter += n_iterations_without_report;
			// iteration with computing L2 change of velocity and stress
			Uchange.save_field();
			const Float Tres = AL.iterate_report_stress_change();
			const Float Ures = Uchange.calculate_field_change();
			residuals_monitor.add_point(++niter,{Ures,Tres});

			if( time_to_print_header.alarm_ringing() ){
				program_output.print("\niteration","|Un+1-Un|_L2","|Tn+1-Tn|_L2");
				program_output.fill_horizontal('-',1.1);
			}
			program_output.print(niter,Ures,Tres);
		} while( (niter<max_iteration) && !residuals_monitor.is_converged() );

		print_solution_convergence_message( residuals_monitor.is_converged() );
		AL.write_results();
		residuals_monitor.save_to_file();
	}

	field adapt_criteria() const
	{return AL.adapt_criteria();}

private:
	VelocityRHSManipulator rhs_manipulator;
	AugmentedLagrangian_basic<VelocityMinimizationSolver> AL;

	int const max_iteration;
	int const n_iterations_without_report;
	L2norm_calculator Uchange;
	ConvergenceMonitor residuals_monitor;
	RecuringAlarm time_to_print_header;

};

#endif /* STANDARDAUGMENTEDLAGRANGIAN_H_ */
