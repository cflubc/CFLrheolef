/*
 * AdaptationLoop.h
 *
 *  Created on: 2013-04-14
 *      Author: ali
 */

#ifndef ADAPTATIONLOOP_H_
#define ADAPTATIONLOOP_H_

#include <cstdlib>
#include <string>
#include <sstream>

#include "rheolef.h"
#include "rheolef/adapt.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "OperatingSystem.h"
#include "adaptationStrategy.h"
#include "PrintArguments.h"



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
	std::string base_name;

public:
	typedef FieldsPool_ FieldsPool;

	AdaptationLoop( XMLConfigFile const& cf,
					FieldsPool& fields,
					DirichletBC& BC
			         ):
		conf(cf),
		strategy( cf.child("Adaptation") ),
		bc(BC),
		base_name(fields.get_geo().name())
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
			OS::chdir_up();
			rheolef::geo const omega = rheolef::adapt( criteria, strategy.opt() );
			FieldsPool fields(conf.child(CFL_FieldsPool_Module),omega,bc);
			Application app(conf,fields,bc);
			CFL_mkresult_folder_and_cd_to_it(i);
			strategy.run_app(app);
			if( i!=strategy.n_adapt )
				criteria = app.adapt_criteria();
		}

		// move files from top folder to corresponding resultx (x adapt number) folder
		OS::chdir_up();
		for(int i=1; i<=strategy.n_adapt; ++i){
			std::stringstream ss;
			print_args(ss,"mv ",base_name,"-",i,"* ",CFL_SaveFolder_BaseName,i);
			system( ss.str().c_str() );
		}
		// move files for first adaptation
		std::stringstream ss;
		print_args(ss,"mv ",base_name,"-* ",CFL_SaveFolder_BaseName,"0");
		system( ss.str().c_str() );

		ss.str("");
		print_args(ss,"mv ",base_name,".* ",CFL_SaveFolder_BaseName,"0");
		system( ss.str().c_str() );
	}
};


#endif /* ADAPTATIONLOOP_H_ */
