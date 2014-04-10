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
#include "OperatingSystem.h"



typedef Problem_MacroBubbleFlowOnset  Problem;
typedef Problem::BC DirichletBoundaryConditions;
typedef Problem::FieldsPool FieldsPool;
typedef Problem::Mesh  Mesh;
typedef NonAdaptiveDriver<Problem::Application,Problem::FieldsPool> Driver;


int main(int argc, char** argv )
{
	using std::string;
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
	string const base_name = meshxml("name");
	if( generate_mesh )
	{
		Mesh create_meshcad( meshxml, base_name );
		string const args = meshxml.get_if_path_exist("command_line_args",string(" "));
		make_geo_from_bamgcad_and_dmn_file( meshxml("npoint",size_t()), base_name, args );
	}
	string const Noplot = "no";
	string const plot_mesh_args = meshxml.get_if_path_exist("plot_mesh_args",Noplot);
	if( plot_mesh_args!=Noplot )
		OS::run_command( "geo -noverbose "+plot_mesh_args+" "+base_name );
	if( exit_after_mesh_gen )
		return EXIT_SUCCESS;

	DirichletBoundaryConditions BC( conf.child(CFL_FieldsPool_Module) );
	Driver::run(conf, base_name, BC);
	timer.stop();

	CFL_print_time_memory_useage(rheolef::dout,timer.get_time_passed());
	std::ofstream  o("useage.info");
	CFL_print_time_memory_useage(o,timer.get_time_passed());
	o.close();

	return EXIT_SUCCESS;
}
