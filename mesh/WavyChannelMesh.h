/*
 * WavyChannelMesh.h
 *
 *  Created on: 2013-06-08
 *      Author: ali
 */

#ifndef WAVYCHANNELMESH_H_
#define WAVYCHANNELMESH_H_

#include <string>
#include <vector>
#include <cstddef>
#include <stdexcept>

#include "ConfigXML.h"
#include "bamgcad.h"
#include "ParametricCurves.h"
#include "ParametricCurveMeshGen.h"


class WavyChannelMesh
{
	typedef std::size_t size_t;

public:

	WavyChannelMesh( XMLConfigFile const& conf, std::string const& base_name ):
		wall( .5*conf("H",double()), conf("L",double()) ),
		paramesh( conf({"ParametricCurve_mesh","max_dTheta"},double()),
				  conf({"ParametricCurve_mesh","max_ds"},double()) )
	{
		std::string const type( conf("type") );
		if( type=="symxy" )
			symxy(base_name);
		else
			throw std::logic_error("Wrong type for ChannelMesh provided");
	}


private:

	wavy_wall wall;
	ParametricCurveMeshGen<double,false,false,true> paramesh;
	std::vector<double> X;
	std::vector<double> Y;


	void symxy( std::string const& base_name )
	{
		paramesh.gen_mesh(wall,-wall.half_of_wavelength,0.,&X,&Y);

		size_t const nwall = X.size();
		size_t const ntrivial_vertices = 4;
		size_t const nvertices = ntrivial_vertices + nwall;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 ",wall.y(0.)," 1\n");
		bamg.print("0 0 2\n");
		bamg.print(-wall.half_of_wavelength," 0  2\n");
		bamg.print(-wall.half_of_wavelength," 1. 1\n");
		for(auto x=X.cbegin(), y=Y.cbegin(); x!=X.cend(); ++x,++y)
			bamg.print(*x," ",*y," 1\n");

		bamg.print_edges_header();
		bamg.print( "1 2 101\n"
				    "2 3 102\n"
				    "3 4 103\n" );
		bamg.print_ordered_edges(ntrivial_vertices,nwall," 104\n");
		bamg.print(nvertices," 1 104\n");
		bamg.close_file();

		std::ofstream fdmn( domain_filename(base_name).c_str() );
		fdmn << "EdgeDomainNames\n"
				"4\n"
				"right\n"
				"bottom\n"
				"left\n"
				"top";
		fdmn.close();
	}
};

#endif /* WAVYCHANNELMESH_H_ */
