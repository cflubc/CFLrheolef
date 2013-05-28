/*
 * adaptationStrategy.h
 *
 *  Created on: 2013-04-16
 *      Author: ali
 */

#ifndef ADAPTATIONSTRATEGY_H_
#define ADAPTATIONSTRATEGY_H_

#include "rheolef.h"

struct FixedStrategy
{
	rheolef::adapt_option_type opts;
	int n_adapt;

	FixedStrategy( const XMLConfigFile& conf )
	{
		conf("max_vertices",&opts.n_vertices_max);
		conf("hmin",&opts.hmin);
		conf("hmax",&opts.hmax);
		conf("hcoef",&opts.hcoef);
		conf("cycles",&n_adapt);
		conf("err",&opts.err);
		conf("additional",&opts.additional);
	}

	template< typename Application >
	void run_app( Application& app )
	{app.run();}

	void set_cycle( int i )
	{}

	const rheolef::adapt_option_type& opt()
	{ return opts; }
};



#endif /* ADAPTATIONSTRATEGY_H_ */
