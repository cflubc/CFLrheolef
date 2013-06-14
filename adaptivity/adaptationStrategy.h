/*
 * adaptationStrategy.h
 *
 *  Created on: 2013-04-16
 *      Author: ali
 */

#ifndef ADAPTATIONSTRATEGY_H_
#define ADAPTATIONSTRATEGY_H_

#include <cstddef>
#include <vector>
#include <string>

#include "rheolef.h"

class FixedStrategy
{
	typedef std::size_t size_t;

public:

	int const n_adapt;

	FixedStrategy( XMLConfigFile const& conf ):
		XML_INIT_VAR(conf,n_adapt,"cycles"),
		hmin(conf,"hmin"),
		hmax(conf,"hmax"),
		hcoef(conf,"hcoef"),
		err(conf,"err"),
		n_vertices_max(conf,"n_vertices_max"),
		ratio(conf,"ratio")
//		errg(conf,"errg")
	{
		opts.additional = conf("additional");
	}

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
		SET_OPT(ratio)
//		SET_OPT(errg)
	#undef SET_OPT
		return opts;
	}


private:

	struct adaptive_parameter {
		typedef std::vector<rheolef::Float> vec;
		vec const p;

		adaptive_parameter() {}

		adaptive_parameter( XMLConfigFile const& conf, cstr const name ):
			p( conf(name,vec()) )
		{}

		rheolef::Float operator[]( size_t const icycle ) const {
			if( icycle<p.size() )
				return p[icycle];
			return p.back();
		}
	};


	rheolef::adapt_option_type opts;
	adaptive_parameter const hmin;
	adaptive_parameter const hmax;
	adaptive_parameter const hcoef;
	adaptive_parameter const err;
	adaptive_parameter const n_vertices_max;
	adaptive_parameter const ratio;
//	adaptive_parameter const errg;
};



#endif /* ADAPTATIONSTRATEGY_H_ */
