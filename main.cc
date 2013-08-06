/*
 * stokes.cc
 *
 *  Created on: 2013-04-07
 *      Author: ali
 */

#include <cstdlib>
#include <stdexcept>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "ApplicationDriver.h"
#include "adaptationCriterions.h"
#include "TimeGauge.h"
#include "Problems.h"


// saving for restart: is a solver adapter? or new wrapper object?
/*
 * restart class: it runs the whole system first, if it is a case
 * of restart, reads the proper geometry and pass to adapt loop
 */

typedef Problem_BubbleEncapsulation   Problem;
typedef Problem::BC  DirichletBoundaryConditions;
typedef Problem::FieldsPool  FieldsPool;
typedef Problem::Mesh  Mesh;
typedef AdaptiveDriver<Problem::Application,Problem::FieldsPool> Driver;


/**/
int main(int argc, char** argv )
{
	rheolef::environment env(argc,argv);

	bool generate_mesh = true;
	bool exit_after_mesh_gen = false;
	for(int i=1; i<argc; ++i)
	{
		cstr arg = argv[i];
		if( cstrcmp(arg,"-no-mesh") )
			generate_mesh = false;
		else if( cstrcmp(arg,"-exit") )
			exit_after_mesh_gen = true;
	}

	XMLConfigFile const conf( argv[1] );
	if( !cstrcmp(Problem::Name,conf("problem")) )
		throw std::logic_error("Name of the problem in xml file doesn't match with source code");

	TimeGauge timer;
	timer.start();

	XMLConfigFile const meshxml = conf.child("Mesh");
	std::string const base_name = meshxml("name");
	if( generate_mesh )
	{
		Mesh create_meshcad( meshxml, base_name );
		make_geo_from_bamgcad_and_dmn_file( base_name,meshxml("command_line_args") );
	}
	plot_mesh(meshxml, base_name);
	if( exit_after_mesh_gen )
		exit(EXIT_SUCCESS);

	DirichletBoundaryConditions BC( conf.child(CFL_FieldsPool_Module) );
	Driver::run(conf, base_name, BC);
	timer.stop();

	CFL_print_time_memory_useage(rheolef::dout,timer.get_time_passed());
	std::ofstream o("useage.info");
	CFL_print_time_memory_useage(o,timer.get_time_passed());
	o.close();

	return EXIT_SUCCESS;
}
/**/

template< typename Application >
struct NonAdaptiveApplication
{
	template< typename Fields, typename BoundaryCondition >
	void setup_and_run( const XMLConfigFile& conf, Fields& fields, BoundaryCondition& BC ){
		Application app(conf,fields);
		app.run();
	}
};

//int main( ){
//	rheolef::idiststream o("wav-1-crit","field");
//	rheolef::field w1, w;
//	o >> w1;
//	o.close();
//
//	o.open("wav-crit","field");
//	o >> w;
//	o.close();
//
//	rheolef::field wi = interpolate(w.get_space(),w1);
//	rheolef::form m(w.get_space(),w.get_space(),"mass");
//	std::cout << "this is the norm " << m(wi-w,wi-w) << '\n';
//}
