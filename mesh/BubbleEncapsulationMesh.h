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
#include <limits>

#include "CFL.h"
#include "ConfigXML.h"
#include "MathUtility.h"
#include "PrintArguments.h"
#include "bamgcad.h"
#include "ParametricCurves.h"
#include "CFLParametricMeshGen.h"
#include "RheolefDomainFile.h"


class BubbleEncapsulationMesh
{
	typedef std::size_t size_t;
	typedef std::string string;
	typedef std::vector<double> vec;

public:

	BubbleEncapsulationMesh( XMLConfigFile const& conf, string const& base_name ):
		type( type_str(conf)[0] ),
		h_wallMesh( conf.get_if_path_exist("hwall",double(0.)) ),
		use_fineMesh_on_wall( 0<h_wallMesh ),
		gen_droplet_mesh( type_str(conf)[1]==string("droplet")  ),
		rx( .5*conf("bubble_length",rx) ),
		ry( .5*conf("bubble_width",ry) ),
		 D( .5*conf("channel_width",D) ),
		 L( .5*conf("channel_length",L) ),
		 A( type=="symxy" ? L*D : 2.*L*D ),
		shape(rx,ry),
		curve_integrator( conf.child("ParametricCurve_mesh") )
	{
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
		size_t  ncorner_vertices = 5;
		if(gen_droplet_mesh)
			++ncorner_vertices;
		size_t const nvertices = nbubble+nwall+ncorner_vertices;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 ",  ry, " 1\n");
		bamg.print("0 ",  D,  " 2\n");
		bamg.print_points(Xwal,Ywal," 2\n");
		bamg.print(-L," ",D,  " 2\n");
		bamg.print(-L,  " 0", " 4\n");
		bamg.print(-rx, " 0", " 4\n");
		bamg.print_points(Xbub,Ybub," 5\n");
		if(gen_droplet_mesh)
			bamg.print("0 0 1\n");


		size_t nedges = nvertices;
		if( gen_droplet_mesh )
			++nedges;
		bamg.print_edges_header(nedges);
		size_t ipoint = 1;
		ipoint = bamg.print_edge(ipoint," 1\n");
		ipoint = bamg.print_ordered_edges(ipoint,nwall+1," 2\n");
		ipoint = bamg.print_edge(ipoint," 3\n");
		ipoint = bamg.print_edge(ipoint," 4\n");
		size_t const bub_head_point_id = ipoint;
		ipoint = bamg.print_ordered_edges(ipoint,nbubble," 5\n");
		bamg.print(ipoint," 1 5\n");

		if( gen_droplet_mesh )
		{
			size_t const last_bub_edge = ipoint;
			++ipoint;
			bamg.print(ipoint," 1 1\n");
			bamg.print(bub_head_point_id," ",ipoint," 4\n");

			bamg.print_subdomain_header(2);
			bamg.print("2 ",last_bub_edge," -1 201\n");
			bamg.print("2 ",last_bub_edge,"  1 202\n");
		}
		bamg.close_file();

		RheolefDomainFile fdmn(base_name);
		if( gen_droplet_mesh ){
			fdmn.print_edge_domains({"right","top","left","bottom","droplet_boundary"});
			fdmn.print_region_domains({"droplet","fluid"});
		}
		else
			fdmn.print_edge_domains({"right","top","left","bottom","bubble"});
		fdmn.close_file();
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
		bamg.print("0 "  , ry ," 1\n");
		bamg.print("0 "  , D  ," 2\n");
		bamg.print_points(Xwal,Ywal," 2\n");
		bamg.print(-L," ", D  ," 2\n");

		for(auto& y:Ywal)
			y*=-1;
		// Ywal=-D for all, no reverse required
		std::reverse( begin(Xwal), end(Xwal) );
		bamg.print(-L," ",-D   ," 4\n");
		bamg.print_points(Xwal,Ywal," 4\n");
		bamg.print("0 "  ,-D   ," 4\n");
		bamg.print("0 "  ,-ry  ," 5\n");
		bamg.print_points(Xbub,Ybub," 6\n");

		size_t nedges = nvertices;
		if( gen_droplet_mesh )
			++nedges;
		bamg.print_edges_header(nedges);
		bamg.print("1 2 1\n");
		size_t iedge = bamg.print_ordered_edges(2,nwall+1," 2\n");
		iedge = bamg.print_edge(iedge," 3\n");
		iedge = bamg.print_ordered_edges(iedge,nwall+1," 4\n");
		iedge = bamg.print_edge(iedge," 1\n");
		iedge = bamg.print_ordered_edges(iedge,nbubble," 5\n");
		bamg.print(iedge," 1 5\n");

		if( gen_droplet_mesh )
		{
			size_t const bub_bottom_point_id = ncorner_vertices+2*nwall;
			bamg.print(bub_bottom_point_id," 1 1\n");
			size_t const& direction_edge_id = iedge;
			bamg.print_subdomain_header(2);
			bamg.print("2 ",direction_edge_id," -1 201\n");
			bamg.print("2 ",direction_edge_id,"  1 202\n");
		}
		bamg.close_file();

		RheolefDomainFile fdmn(base_name);
		if( gen_droplet_mesh ){
			fdmn.print_edge_domains({"right","top","left","bottom","droplet_boundary"});
			fdmn.print_region_domains({"droplet","fluid"});
		}
		else
			fdmn.print_edge_domains({"right","top","left","bottom","bubble"});
		fdmn.close_file();
	}

	double area() const
	{return A;}

private:

	std::vector<string> type_str( XMLConfigFile const& conf ){
		std::vector<string> v;
		conf("type",&v);
		return v;
	}

	std::string const type;
	double const h_wallMesh;
	bool const use_fineMesh_on_wall;
	bool const gen_droplet_mesh;

	double const rx;
	double const ry;
	double const  D;
	double const  L;
	double const A;

	shape_ellipse shape;
	CFLCurveIntegrator curve_integrator;
};


#endif /* BUBBLEENCAPSULATIONMESH_H_ */
