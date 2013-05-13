/*
 * StandardAugmentedLagrangian.h
 *
 *  Created on: 2013-05-02
 *      Author: ali
 */

#ifndef STANDARDAUGMENTEDLAGRANGIAN_H_
#define STANDARDAUGMENTEDLAGRANGIAN_H_

#include "rheolef.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "OutputFormatting.h"
#include "ConvergenceMonitor.h"
#include "AugmentedLagrangian_basic.h"


template< typename FlowSolver >
class StandardAugmentedLagrangian
{
	typedef rheolef::Float Float;
	typedef rheolef::field field;

public:
	template< typename FieldsPool, typename DirichletBC >
	StandardAugmentedLagrangian( const XMLConfigFile& conf,
								 FieldsPool& fields,      //const rheolef::geo& omega,
								 DirichletBC BC
								 ):
		AL(conf,fields,BC),
		TminusaG_rhs(fields.Xh, 0.),
		n_report( conf.atoi_if_exist("reports_frequency",10) ),
		max_iteration( conf.atoi("max_iteration") ),
		time_to_print_header( conf.atoi_if_exist("report_header_reprint_frequency",30) ),
		residuals_monitor("UTconvergence", conf.atof("convergence_limit"), {"|Un+1-Un|","|Tn+1-Tn|"}),
		Uchange(fields.Uh)
	{}

	void run()
	{
		auto program_output( make_column_output(std::cout,10,16,16) );
		int niter(0);
//		L2norm_calculator Uchange( AL.Uh() );

		do {
			// iterations without reporting
			for(int i=0; i<n_report-1; ++i){
				AL.solve(TminusaG_rhs);
				AL.contribute_to_rhs_fast(TminusaG_rhs);
				++niter;
			}
			// iteration with computing L2 change of velocity and stress
			Uchange.save_field();
			AL.solve(TminusaG_rhs);
			const Float Tres = AL.contribute_to_rhs_report_stress_change(TminusaG_rhs);
			const Float Ures = Uchange.calculate_field_change();
			++niter;
			residuals_monitor.add_point(niter,{Ures,Tres});

			if( time_to_print_header.alarm_ringing() ){
				program_output.print("\niteration","|Un+1-Un|_L2","|Tn+1-Tn|_L2");
				program_output.fill_horizontal('-',1.1);
			}
			program_output.print(niter,Ures,Tres);
		} while( !residuals_monitor.is_converged() && (niter<max_iteration) );

		if( niter<max_iteration )
			std::cout << "\nThe solution converged... :-)\n";
		else
			std::cout << "\nMax limit of iterations reached, stopping...\n";

		AL.write_results();
		residuals_monitor.save_to_file();
	}

	field adapt_criteria() const
	{ return AL.adapt_criteria(); }

private:
	AugmentedLagrangian_basic<FlowSolver> AL;
	rheolef::field TminusaG_rhs;

	int n_report;
	int max_iteration;
	RecuringAlarm time_to_print_header;
	ConvergenceMonitor residuals_monitor;
	L2norm_calculator Uchange;

};

#endif /* STANDARDAUGMENTEDLAGRANGIAN_H_ */
