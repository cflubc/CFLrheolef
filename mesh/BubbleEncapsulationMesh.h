/*
 * BubbleEncapsulationMesh.h
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */

#ifndef BUBBLEENCAPSULATIONMESH_H_
#define BUBBLEENCAPSULATIONMESH_H_


#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdexcept>

#include "CFL.h"
#include "ConfigXML.h"
#include "PrintArguments.h"
#include "bamgcad.h"
#include "ParametricCurves.h"
#include "ParametricCurveMeshGen.h"


class BubbleEncapsulationMesh
{
	typedef std::string string;

public:
	BubbleEncapsulationMesh( XMLConfigFile const& conf, string const& base_name ):
		rx( .5*conf.atof("bubble_length") ),
		ry( .5*conf.atof("bubble_width") ),
		 D( .5*conf.atof("channel_width") ),
		 L( .5*conf.atof("channel_length") ),
		shape(rx,ry),
		paramesh( conf.atof("bubble_curve_dteta_limit") )
	{
		string const type( conf("type") );
		if( type=="symx" )
			symx(base_name);
		else if( type=="symxy" )
			symxy(base_name);
		else
			throw std::logic_error("Wrong type of Bubble mesh");
	}

	void symxy( string const& base_name )
	{
		paramesh.gen_mesh(shape,1.5*PI,2.*PI,&X,&Y);

		size_t const nbubble = X.size();
		size_t const ntrivial_vertices = 5;
		size_t const nvertices = nbubble+ntrivial_vertices;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 ",  ry, " 1\n");
		bamg.print("0 ",  D,  " 2\n");
		bamg.print(-L," ",D,  " 2\n");
		bamg.print(-L,  " 0", " 4\n");
		bamg.print(-rx, " 0", " 4\n");
		for(auto x=X.begin(), y=Y.begin(); x!=X.end(); ++x,++y )
			bamg.print(*x," ",*y," 5\n");

		bamg.print_edges_header();
		bamg.print( "1 2 1\n"
				    "2 3 2\n"
				    "3 4 3\n"
				    "4 5 4\n" );
		bamg.print_ordered_edges(5,nbubble," 5\n");
		bamg.print(nbubble+5," 1 5\n");
		bamg.close_file();

		std::ofstream fdmn( domain_filename(base_name).c_str() );
		fdmn << "EdgeDomainNames\n"
				"5\n"
				"right_top\n"
				"top\n"
				"left\n"
				"bottom\n"
				"bubble\n";
		fdmn.close();
	}

	void symx( string const& base_name )
	{
		paramesh.gen_mesh(shape,PI,2.*PI,&X,&Y);

		size_t const nbubble = X.size();
		size_t const ntrivial_vertices = 6;
		size_t const nvertices = nbubble+ntrivial_vertices;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 "  , ry  ," 1\n");
		bamg.print("0 "  , D   ," 2\n");
		bamg.print(-L," ", D   ," 2\n");
		bamg.print(-L," ",-D   ," 4\n");
		bamg.print("0 "  ,-D   ," 4\n");
		bamg.print("0 "  ,-ry  ," 5\n");
		for(auto x=X.begin(), y=Y.begin(); x!=X.end(); ++x,++y)
			bamg.print(*x," ",*y," 6\n");

		bamg.print_edges_header();
		bamg.print(
				   "1 2 1\n"
				   "2 3 2\n"
				   "3 4 3\n"
				   "4 5 4\n"
				   "5 6 5\n" );
		bamg.print_ordered_edges(6,nbubble," 6\n");
		bamg.print(nbubble+6," 1 6\n");
		bamg.close_file();

		std::ofstream fdmn( domain_filename(base_name).c_str() );
		fdmn << "EdgeDomainNames\n"
				"6\n"
				"right_top\n"
				"top\n"
				"left\n"
				"bottom\n"
				"right_bottom\n"
				"bubble\n";
		fdmn.close();
	}

private:
	double const rx;
	double const ry;
	double const  D;
	double const  L;

	shape_ellipse shape;
	ParametricCurveMeshGen<double,ParametricMeshOpts::exclude_ends> paramesh;
	std::vector<double> X;
	std::vector<double> Y;

};


#endif /* BUBBLEENCAPSULATIONMESH_H_ */
