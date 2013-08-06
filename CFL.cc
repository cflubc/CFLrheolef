/*
 * CFL.cc
 *
 *  Created on: 2013-04-17
 *      Author: ali
 */


#include <cstring>
#include <cassert>
#include <cstdio>

#include <string>
#include <stdexcept>

#include "CFL.h"
#include "OperatingSystem.h"
#include "ConfigXML.h"


std::string derivative_approx( const std::string& approx )
{
	if( approx=="P2" )
		return "P1d";
	else if( approx=="P1" )
		return "P0";
	else
		throw std::logic_error("Derivative of "+approx+" FEM elements is not implemented");
}


void assert_equal( const rheolef::field& f1, const rheolef::field& f2 )
{
	assert( f1.ndof()==f2.ndof() );
	assert( f1.size()==f2.size() );
}

void print_solution_convergence_message( bool converged )
{
	if( converged )
		printf("\nThe solution converged... :-)\n");
	else
		printf("\nMax limit of iterations reached, stopping...\n");
}

void CFL_mkresult_folder_and_cd_to_it( int iadapt )
{
	std::string const name(CFL_SaveFolder_BaseName+std::to_string(iadapt));
	OS::mkdir(name);
	OS::changedir(name);
}


void plot_mesh( XMLConfigFile const& conf, std::string const& base_name )
{
	cstr const no_disp = "no";
	cstr const args = conf.return_txt_if_exist({"plot_mesh_args"},no_disp);
	if( strcmp(no_disp,args)==0 )
		return;
	OS::run_command( "geo -noverbose "+base_name+" "+args );
}

