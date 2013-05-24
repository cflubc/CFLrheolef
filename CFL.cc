/*
 * CFL.cc
 *
 *  Created on: 2013-04-17
 *      Author: ali
 */


#include <string>
#include <cassert>
#include <stdexcept>
#include <cstdio>

#include "CFL.h"
#include "OperatingSystem.h"


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
		printf("\nThe solution converged... :-)\n\n");
	else
		printf("\nMax limit of iterations reached, stopping...\n\n");
}

void CFL_mkresult_folder_and_cd_to_it( int iadapt )
{
	std::string const name("result"+std::to_string(iadapt));
	OS::mkdir(name);
	OS::changedir(name);
}

