/*
 * GenericIteration.h
 *
 *  Created on: Sep 12, 2013
 *      Author: ali
 */

#ifndef GENERICITERATION_H_
#define GENERICITERATION_H_

#include <cstddef>
#include "rheolef/compiler.h"
#include "ConfigXML.h"


class GenericIteration
{
	typedef rheolef::Float Float;
	typedef std::size_t size_t;

public:

	GenericIteration( const XMLConfigFile& conf ):
		XML_INIT_VAR(conf,convergence_limit,"convergence_limit"),
		XML_INIT_VAR(conf,max_iteration,"max_iteration"),
		XML_INIT_VAR(conf,min_iteration,"min_iteration"),
		XML_INIT_VAR(conf,report_frequency,"reports_frequency")
	{}


	template< typename iteration >
	void operator()( iteration& ITR )
	{
		size_t niter = 0;
		bool is_converged = false;
		do
		{
			Float residual;
			niter += report_frequency;
			for(size_t i=0; i<report_frequency-1; ++i)
				ITR.iterate();
			ITR.iterate_report(niter,residual);
			is_converged = residual<convergence_limit;
		} while(
				(niter<min_iteration) || (niter<max_iteration && !is_converged)
			   );
		print_solution_convergence_message(is_converged);
	}

	Float  const convergence_limit;
	size_t const max_iteration;
	size_t const min_iteration;
	size_t const report_frequency;
};



#endif /* GENERICITERATION_H_ */
