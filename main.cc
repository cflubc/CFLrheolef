/*
 * stokes.cc
 *
 *  Created on: 2013-04-07
 *      Author: ali
 */

#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "adaptationLoop.h"
#include "adaptationCriterions.h"
#include "TimeGauge.h"
#include "Problems.h"


// saving for restart: is a solver adapter? or new wrapper object?
/*
 * restart class: it runs the whole system first, if it is a case
 * of restart, reads the proper geometry and pass to adapt loop
 */

typedef Problem_AugLag_BubbleEncapsulation  Problem;
typedef Problem::BC  DirichletBoundaryConditions;
typedef Problem::FieldsPool  FieldsPool;
typedef Problem::Mesh  Mesh;
//typedef Problem::Application  Application;
typedef AdaptationLoop<Problem::Application,Problem::FieldsPool,Problem::BC>  Application;



int main(int argc, char** argv )
{
	rheolef::environment env(argc,argv);

	XMLConfigFile const conf( argv[1] );
	if( strcmp(Problem::Name,conf("problem"))!=0 )
		throw std::logic_error("Name of the problem in xml file doesn't match with source code");

	TimeGauge timer;
	timer.start();

	XMLConfigFile const meshxml = conf.child("mesh");
	std::string const base_name = meshxml("name");
	std::string const geo_file = geo_filename(base_name);
	Mesh create_meshcad( meshxml, base_name );
	make_geo_from_bamgcad_and_dmn_file( base_name,meshxml("command_line_args") );
	plot_mesh(meshxml, geo_file);
	if(argc==3)
		exit(0);

	rheolef::geo omega(geo_file);
	XMLConfigFile const FE( conf.child(CFL_FieldsPool_Module) );
	DirichletBoundaryConditions BC(FE);
	FieldsPool fields(FE,omega,BC);

	CFL_mkresult_folder_and_cd_to_it(0);
	Application app(conf,fields,BC);
	app.run();

	timer.stop();

	CFL_print_time_memory_useage(rheolef::dout,timer.get_time_passed());
	std::ofstream o("useage.info");
	CFL_print_time_memory_useage(o,timer.get_time_passed());
	o.close();

	return EXIT_SUCCESS;
}


template< typename Application >
struct NonAdaptiveApplication
{
	template< typename Fields, typename BoundaryCondition >
	void setup_and_run( const XMLConfigFile& conf, Fields& fields, BoundaryCondition& BC ){
		Application app(conf,fields);
		app.run();
	}
};
