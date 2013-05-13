/*
 * AdaptationLoop.h
 *
 *  Created on: 2013-04-14
 *      Author: ali
 */

#ifndef ADAPTATIONLOOP_H_
#define ADAPTATIONLOOP_H_

#include <cassert>
#include "rheolef.h"
#include "rheolef/adapt.h"
#include "rheolef/diststream.h"

#include "ConfigXML.h"
#include "adaptationStrategy.h"



template< typename Application,
	      typename FieldsPool,
	      typename DirichletBC,
    	  typename AdaptStrategy = FixedStrategy >
class AdaptationLoop
{
	const XMLConfigFile& conf;
	AdaptStrategy strategy;
	DirichletBC& bc;
	rheolef::field criteria;

public:
	AdaptationLoop( const XMLConfigFile& cf,
					FieldsPool fields,
					DirichletBC& BC
			         ):
		conf(cf),
		strategy( cf.child("Adaptation") ),
		bc(BC)
	{
		strategy.set_cycle(0);
		Application app(conf,fields,bc);
		strategy.run_app(app);
		criteria = app.adapt_criteria();
	}


	void run()
	{
		for(int i=1; i<=strategy.n_adapt; ++i)
		{
			strategy.set_cycle(i);
			const rheolef::geo& omega = rheolef::adapt( criteria, strategy.opt() );
			FieldsPool fields(conf,omega,bc);
			Application app(conf,fields,bc);
			strategy.run_app(app);
			if( i!=strategy.n_adapt )
				criteria = app.adapt_criteria();
		}
	}
};


#endif /* ADAPTATIONLOOP_H_ */
