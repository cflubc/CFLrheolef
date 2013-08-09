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
#include "MathUtility.h"
#include "bamgcad.h"
#include "ParametricCurves.h"
#include "CFLParametricMeshGen.h"


class WavyChannelMesh
{
	typedef std::size_t size_t;

public:

	WavyChannelMesh( XMLConfigFile const& conf, std::string const& base_name ):
		inlet_length( conf.get_if_path_exist("Inlet",0) ),
		wall( .5*conf("H",double()), conf("L",double()) ),
		curve_integrator( conf.child("ParametricCurve_mesh") )
	{
		if( inlet_length<0 )
			throw std::logic_error("Inlet length of channel should be positive!");
		X.reserve(300);
		Y.reserve(300);
		std::string const type( conf("type") );
		if( type=="symxy" )
			symxy(base_name);
		else
			throw std::logic_error("Wrong type for ChannelMesh provided");
	}


private:

	void symxy( std::string const& base_name )
	{
		ExclusiveInterval<double> const range(-wall.half_of_wavelength,0.);
		gen_parametric_curve_mesh(wall,curve_integrator,range,&X,&Y);

		bool const has_inlet = (0<inlet_length);
		size_t const nwall = X.size();
		size_t const ntrivial_vertices = has_inlet ? 5:4;
		size_t const nvertices = ntrivial_vertices + nwall;

		bamgcad bamg( nvertices, base_name );
		bamg.print("0 ",wall.y(0.)," 1\n");
		bamg.print("0 0 2\n");
		double const xinlet = -( wall.half_of_wavelength + inlet_length );
		bamg.print(xinlet," 0 2\n");
		bamg.print(xinlet," 1 1\n");
		if( has_inlet )
			bamg.print(-wall.half_of_wavelength," 1 1\n");
		for(auto x=X.cbegin(), y=Y.cbegin(); x!=X.cend(); ++x,++y)
			bamg.print(*x," ",*y," 1\n");

		bamg.print_edges_header();
		bamg.print( "1 2 101\n"
				    "2 3 102\n"
				    "3 4 103\n" );
		size_t const npoint_on_top = has_inlet ? nwall+1:nwall;
		bamg.print_ordered_edges(4,npoint_on_top," 104\n");
		bamg.print(nvertices," 1 104\n");
		bamg.close_file();

		RheolefDomainFile fdmn(base_name);
		fdmn.print_edge_domains({"right","bottom","left","top"});
		fdmn.close_file();
	}

	double const inlet_length;
	wavy_wall wall;
	CFLCurveIntegrator curve_integrator;
	std::vector<double> X;
	std::vector<double> Y;
};

#endif /* WAVYCHANNELMESH_H_ */
