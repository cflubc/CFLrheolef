/*
 * bamgcad.cc
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */


#include <cstdlib>
#include <iostream>
#include <sstream>

#include "CFL.h"
#include "bamgcad.h"

void make_geo_from_bamgcad_and_dmn_file(
								std::string const& base_name,
								std::string const& other_commandline_args )
{
	std::string const bamgfile = bamgmesh_filename(base_name);
	std::ostringstream os;
	print_args(os,
		"bamg -g ",bamgcad_filename(base_name)," -o ",bamgfile," ",other_commandline_args);
	println_args( std::cout,"[Exec] ",os.str() );
	system( os.str().c_str() );


	std::ostringstream os2;
	print_args(os2,
			"bamg2geo ",domain_filename(base_name)," ",bamgfile," > ",geo_filename(base_name));
	println_args( std::cout, "[Exec] ",os2.str() );
	system( os2.str().c_str() );
}


