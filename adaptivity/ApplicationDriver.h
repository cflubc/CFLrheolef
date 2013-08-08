/*
 * ApplicationDriver.h
 *
 *  Created on: Aug 5, 2013
 *      Author: ali
 */

#ifndef APPLICATIONDRIVER_H_
#define APPLICATIONDRIVER_H_

#include <cstdlib>
#include <string>
#include "rheolef.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "OperatingSystem.h"
#include "adaptationStrategy.h"



template< typename Application, typename FieldsPool >
struct NonAdaptiveDriver
{
	template< typename BoundaryCondition >
	static void run( XMLConfigFile const& conf, std::string const& geo_basename, BoundaryCondition& BC )
	{
		rheolef::geo omega(geo_basename);
		FieldsPool fields(conf.child(CFL_FieldsPool_Module), omega, BC);
		Application app(conf,fields,BC);
		app.run();
	}
};


template< typename Application, typename FieldsPool >
struct AdaptiveDriver
{
	template< typename BoundaryCondition >
	static void run( XMLConfigFile const& conf, std::string const& geo_basename, BoundaryCondition const& BC )
	{
		FixedStrategy strategy( conf.child("Adaptation") );
		std::string const base_name(geo_basename);

		rheolef::geo omega(geo_basename);
		for(int i=0; i<=strategy.n_adapt; ++i)
		{
			FieldsPool fields(conf.child(CFL_FieldsPool_Module),omega,BC);
			Application app(conf,fields,BC);
			CFL_mkresult_folder_and_cd_to_it(i);
			strategy.run_app(app);
			OS::chdir_up();
			if( i!=strategy.n_adapt ){
				rheolef::field const criteria = app.adapt_criteria();
				omega = rheolef::adapt( criteria, strategy.opt(i) );
			}
		}

		// move files from top folder to corresponding adaptaion folder
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


#endif /* APPLICATIONDRIVER_H_ */
