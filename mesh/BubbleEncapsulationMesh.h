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
	typedef std::size_t size_t;
	typedef std::string string;

public:

	BubbleEncapsulationMesh( XMLConfigFile const& conf, string const& base_name ):
		rx( .5*XML_VAL(conf,rx,"bubble_length") ),
		ry( .5*XML_VAL(conf,ry,"bubble_width") ),
		 D( .5*XML_VAL(conf,D,"channel_width") ),
		 L( .5*XML_VAL(conf,L,"channel_length") ),
		shape(rx,ry),
		paramesh( conf({"ParametricCurve_mesh","max_dTheta"},rx),
				  conf({"ParametricCurve_mesh","max_ds"},rx) )
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
//		paramesh.gen_mesh(shape,1.5*PI,2.*PI,&X,&Y);

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
//		paramesh.gen_mesh(shape,PI,2.*PI,&X,&Y);
		double const tbeg = PI;
		double const tend = 2.*PI;
		size_t const Np = 30;
		double const dt = (tend-tbeg)/(Np+1);
		std::vector<double> points(Np);
		X.resize(Np);
		Y.resize(Np);
		double t = tbeg;
		for(size_t i=0; i<Np; ++i){
			t += dt;
			X[i] = shape.x(t);
			Y[i] = shape.y(t);
		}


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
	ParametricCurveMeshGen<double,false,false,true> paramesh;
	std::vector<double> X;
	std::vector<double> Y;

};


#endif /* BUBBLEENCAPSULATIONMESH_H_ */
