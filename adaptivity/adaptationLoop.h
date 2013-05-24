/*
 * AdaptationLoop.h
 *
 *  Created on: 2013-04-14
 *      Author: ali
 */

#ifndef ADAPTATIONLOOP_H_
#define ADAPTATIONLOOP_H_

#include "rheolef.h"
#include "rheolef/adapt.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "OperatingSystem.h"
#include "adaptationStrategy.h"



template< typename Application,
	      typename FieldsPool_,
	      typename DirichletBC,
    	  typename AdaptStrategy = FixedStrategy >
class AdaptationLoop
{
	XMLConfigFile const& conf;
	AdaptStrategy strategy;
	DirichletBC& bc;
	rheolef::field criteria;

public:
	typedef FieldsPool_ FieldsPool;

	AdaptationLoop( XMLConfigFile const& cf,
					FieldsPool& fields,
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
			OS::changedir("..");
			rheolef::geo const omega = rheolef::adapt( criteria, strategy.opt() );
			FieldsPool fields(conf.child(FieldsPool_Module),omega,bc);
			Application app(conf,fields,bc);
			CFL_mkresult_folder_and_cd_to_it(i);
			strategy.run_app(app);
			if( i!=strategy.n_adapt )
				criteria = app.adapt_criteria();
		}
	}
};


#endif /* ADAPTATIONLOOP_H_ */
