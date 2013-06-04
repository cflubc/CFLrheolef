/*
 * adaptationStrategy.h
 *
 *  Created on: 2013-04-16
 *      Author: ali
 */

#ifndef ADAPTATIONSTRATEGY_H_
#define ADAPTATIONSTRATEGY_H_

#include <cstdlib>
#include <vector>
#include <string>

#include "rheolef.h"

class FixedStrategy
{
public:
	int const n_adapt;

	FixedStrategy( XMLConfigFile const& conf ):
		XML_INIT_VAR(conf,n_adapt,"cycles"),
		hmin(conf,"hmin"),
		hmax(conf,"hmax"),
		hcoef(conf,"hcoef"),
		err(conf,"err"),
		n_vertices_max(conf,"n_vertices_max")
	{}

	template< typename Application >
	void run_app( Application& app )
	{app.run();}


	rheolef::adapt_option_type const& opt( size_t const icycle )
	{
	#define SET_OPT(name) opts.name = name[icycle];
		SET_OPT(hmin)
		SET_OPT(hmax)
		SET_OPT(hcoef)
		SET_OPT(err)
		SET_OPT(n_vertices_max)
	#undef SET_OPT
		return opts;
	}

private:
	struct adaptive_parameter {
		typedef std::vector<rheolef::Float> vec;
		vec const p;

		adaptive_parameter( XMLConfigFile const& conf, cstr const name ):
			p( conf(name,vec()) )
		{}

		rheolef::Float operator[]( size_t icycle ) const {
			if( icycle<p.size() )
				return p[icycle];
			return p.back();
		}
	};


	rheolef::adapt_option_type opts;
	adaptive_parameter hmin;
	adaptive_parameter hmax;
	adaptive_parameter hcoef;
	adaptive_parameter err;
	adaptive_parameter n_vertices_max;
};



#endif /* ADAPTATIONSTRATEGY_H_ */
