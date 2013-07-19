/*
 * BubbleEncapsulationMesh.h
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */

#ifndef BUBBLEENCAPSULATIONMESH_H_
#define BUBBLEENCAPSULATIONMESH_H_


#include <cstddef>
#include <cmath>
#include <cstring>
#include <vector>
#include <fstream>
#include <stdexcept>

#include "CFL.h"
#include "ConfigXML.h"
#include "MathUtility.h"
#include "PrintArguments.h"
#include "bamgcad.h"
#include "ParametricCurves.h"
#include "CFLParametricMeshGen.h"


class BubbleEncapsulationMesh
{
	typedef std::size_t size_t;
	typedef std::string string;
	typedef std::vector<double> vec;

public:

	BubbleEncapsulationMesh( XMLConfigFile const& conf, string const& base_name ):
		rx( .5*conf("bubble_length",rx) ),
		ry( .5*conf("bubble_width",ry) ),
		 D( .5*conf("channel_width",D) ),
		 L( .5*conf("channel_length",L) ),
		use_fineMesh_on_wall( conf.get_if_path_exist({"use_fineMesh_on_wall"},string("no"))==string("yes") ),
		h_wallMesh( use_fineMesh_on_wall ? conf("hwall",h_wallMesh) : 1e7 ),
		shape(rx,ry),
		curve_integrator( conf.child("ParametricCurve_mesh") )
	{
		string const type( conf("type") );
		if( type=="symx" )
			symx(base_name);
		else if( type=="symxy" )
			symxy(base_name);
		else
			throw std::logic_error("Wrong type of Bubble mesh, either symx/symxy");
	}

	void symxy( string const& base_name )
	{
		ExclusiveInterval<double> range_bub(1.5*PI,2.*PI);
		vec Xbub, Ybub;
		gen_parametric_curve_mesh(shape,curve_integrator,range_bub,&Xbub,&Ybub);

		straight_line wall(0,D,-L,D);
		vec Xwal, Ywal;
		if( use_fineMesh_on_wall )
			gen_straight_line_mesh(wall,h_wallMesh,&Xwal,&Ywal,range_bub);

		size_t const nbubble = Xbub.size();
		size_t const nwall = Xwal.size();
		size_t const ncorner_vertices = 5;
		size_t const nvertices = nbubble+nwall+ncorner_vertices;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 ",  ry, " 1\n");
		bamg.print("0 ",  D,  " 2\n");
		bamg.print_points(Xwal,Ywal," 2\n");
		bamg.print(-L," ",D,  " 2\n");
		bamg.print(-L,  " 0", " 4\n");
		bamg.print(-rx, " 0", " 4\n");
		bamg.print_points(Xbub,Ybub," 5\n");

		bamg.print_edges_header();
		bamg.print("1 2 1\n");
		size_t
		iedge = bamg.print_ordered_edges(2,nwall+1," 2\n");
		iedge = bamg.print_edge(iedge," 3\n");
		iedge = bamg.print_edge(iedge," 4\n");
		iedge = bamg.print_ordered_edges(iedge,nbubble," 5\n");
		bamg.print(iedge," 1 5\n");
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
		vec Xbub, Ybub;
		ExclusiveInterval<double> range_bub(PI,2.*PI);
		gen_parametric_curve_mesh(shape,curve_integrator,range_bub,&Xbub,&Ybub);

		straight_line wall(0,D,-L,D);
		vec Xwal, Ywal;
		if( use_fineMesh_on_wall )
			gen_straight_line_mesh(wall,h_wallMesh,&Xwal,&Ywal,range_bub);

		size_t const nbubble = Xbub.size();
		size_t const nwall = Xwal.size();
		size_t const ncorner_vertices = 6;
		size_t const nvertices = nbubble+ncorner_vertices+2*nwall;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 "  , ry  ," 1\n");
		bamg.print("0 "  , D   ," 2\n");
		bamg.print_points(Xwal,Ywal," 2\n");
		bamg.print(-L," ", D   ," 2\n");

		for(auto& y:Ywal)
			y*=-1;
		// Ywal=-D for all, no reverse required
		std::reverse( begin(Xwal), end(Xwal) );
		bamg.print(-L," ",-D   ," 4\n");
		bamg.print_points(Xwal,Ywal," 4\n");
		bamg.print("0 "  ,-D   ," 4\n");
		bamg.print("0 "  ,-ry  ," 5\n");
		bamg.print_points(Xbub,Ybub," 6\n");

		bamg.print_edges_header();
		bamg.print("1 2 1\n");
		size_t iedge = bamg.print_ordered_edges(2,nwall+1," 2\n");
		iedge = bamg.print_edge(iedge," 3\n");
		iedge = bamg.print_ordered_edges(iedge,nwall+1," 4\n");
		iedge = bamg.print_edge(iedge," 5\n");
		iedge = bamg.print_ordered_edges(iedge,nbubble," 6\n");
		bamg.print(iedge," 1 6\n");
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
	bool const use_fineMesh_on_wall;
	double const h_wallMesh;

	shape_ellipse shape;
	CFLCurveIntegrator curve_integrator;
};


#endif /* BUBBLEENCAPSULATIONMESH_H_ */
