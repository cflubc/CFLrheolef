/*
 * bamgcad.cc
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */


#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

#include "CFL.h"
#include "bamgcad.h"
#include "PrintArguments.h"
#include "MathUtility.h"


inline static void
extend_range( double const& x, double& min, double& max )
{
	if( x<min )
		min = x;
	if( max<x )
		max = x;
}

void make_geo_from_bamgcad_and_dmn_file(
								std::size_t const& npoints,
								std::string const& base_name,
								std::string const& other_commandline_args )
{
	using namespace std;
	string const cadfile = bamgcad_filename(base_name);

	// find the bounding box of mesh
	ifstream o(cadfile);
	for(string w; w!="Vertices"; o>>w);
	size_t nvertices;
	o>>nvertices;
	double xmin,xmax, ymin,ymax;
	string tmp;
	// read the first two points
	extract_args(o,xmin,ymin,tmp,xmax,ymax,tmp);
	if( xmax<xmin )
		swap(xmin,xmax);
	if( ymax<ymin )
		swap(ymin,ymax);
	for(size_t i=0; i<nvertices-2; ++i){
		double x,y;
		extract_args(o,x,y,tmp);
		extend_range(x,xmin,xmax);
		extend_range(y,ymin,ymax);
	}
	o.close();
	double const h = sqrt( 3.*(ymax-ymin)*(xmax-xmin)/(npoints*PI) );


	string const bamgfile = bamgmesh_filename(base_name);
	ostringstream os;
	print_args(os,
		"bamg -g ",cadfile," -o ",bamgfile," -hmax ",h," -splitpbedge ",other_commandline_args);
	println_args( cout,"[Exec] ",os.str() );
	system( os.str().c_str() );


	ostringstream os2;
	print_args(os2,
			"bamg2geo ",domain_filename(base_name)," ",bamgfile," > ",geo_filename(base_name));
	print_args( cout, "[Exec] ",os2.str(), "\n\n" );
	system( os2.str().c_str() );
}


