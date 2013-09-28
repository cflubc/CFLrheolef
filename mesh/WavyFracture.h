/*
 * WavyFracture.h
 *
 *  Created on: Aug 20, 2013
 *      Author: steven
 */

#ifndef WAVYFRACTURE_H_
#define WAVYFRACTURE_H_

#include <string>
#include <iostream>
#include <vector>
#include <cstddef>
#include <stdexcept>

#include "ConfigXML.h"
#include "CFLParametricMeshGen.h"
#include "PrintArguments.h"
#include "bamgcad.h"



struct fracture_curve
{
	fracture_curve( XMLConfigFile const& conf ):
		XML_INIT_VAR(conf,cycles,"cycles"),
		XML_INIT_VAR(conf,amplitude,"amplitude"),
		XML_INIT_VAR(conf,phase,"phase"),
		length_half( conf("L",length_half)/2. )

	{
		if (cycles%2==0)
			is_even=true;
		else
			is_even=false;
		if( !(0<=phase && phase<=1.) )
			throw std::logic_error("The phase should be in range [0 1]");
		//add a bunch of other refinements here
	}

	const int cycles;
	const double amplitude;
	const double phase;
	const double length_half;
	bool is_even;

	double scaled_t( double const& t ) const
	{return PI*cycles/length_half*t;}
	double x( const double& t ) const
	{return t;}
	double y(const double & t) const
	{
		if (is_even==false)
			return amplitude*(1.+cos(scaled_t(t)))+0.5;
		else
			return amplitude*(1.-cos(scaled_t(t)))+0.5;
	}
	double dx(const double & t) const
	{return 1;}
	double ddx(const double & t) const
	{return 0;}
	double dy(const double & t) const
	{
		if (is_even==false)
			return amplitude*PI*cycles/length_half*-sin(scaled_t(t));
		else
			return amplitude*PI*cycles/length_half*sin(scaled_t(t));
	}
	double ddy(const double & t) const
	{
		if (is_even==false)
			return -amplitude*pow(PI*cycles/length_half,2)*cos(scaled_t(t));
		else
			return amplitude*pow(PI*cycles/length_half,2)*cos(scaled_t(t));
	}
};



class WavyFracture
{
	typedef std::size_t size_t;
	fracture_curve fracture;
	CFLCurveIntegrator curve_integrator;
	double const entry_length;

public:

	WavyFracture(XMLConfigFile const& conf, std::string const& base_name):
		fracture (conf),
		curve_integrator( conf.child("ParametricCurve_mesh") ),
		XML_INIT_VAR(conf,entry_length,"entry_length")
	{
		std::vector <double> X;
		std::vector <double> Y;
		X.reserve(300);
		Y.reserve(300);

		ExclusiveInterval<double> const range(-fracture.length_half,fracture.length_half);
		gen_parametric_curve_mesh(fracture,curve_integrator,range,&X,&Y);

		size_t const nfrac = X.size();
		size_t const n_corners = 8;
		size_t const ntotal=2*nfrac+n_corners;
		bamgcad bamg(ntotal,base_name);

		double const shift_in_x = fracture.phase*fracture.length_half/fracture.cycles;
		double const left_most_x = fracture.length_half + shift_in_x + entry_length;
		double const right_most_x = fracture.length_half+entry_length;
		bamg.print(-left_most_x," 0.5 1\n");
		bamg.print(-fracture.length_half," 0.5 1\n");

		for(auto x=X.begin(), y=Y.begin(); x!=X.end(); ++x,++y)
			bamg.print(*x," ",*y," 1\n");

		bamg.print(fracture.length_half, " 0.5 1 \n");
		bamg.print(right_most_x, "  0.5 2\n");
		bamg.print(right_most_x, " -0.5 2\n");
		bamg.print(fracture.length_half-shift_in_x," -0.5 1\n");

		auto Xrev = X;
		auto Yrev = Y;
		for (auto yr=begin(Yrev), y=begin(Y), yr_end=end(Yrev); yr!=yr_end; ++yr, ++y)
			*yr=-*y;
		std::reverse(begin(Xrev),end(Xrev));
		std::reverse(begin(Yrev),end(Yrev));

		for(auto x=begin(Xrev), x_end=end(Xrev), y=begin(Yrev); x!=x_end; ++x,++y)
				bamg.print(*x-shift_in_x," ",*y," 1\n");

		bamg.print(-fracture.length_half-shift_in_x, " -0.5 1\n");
		bamg.print(-left_most_x, " -0.5 1\n");


		int current_pos;
		bamg.print_edges_header();
		current_pos=bamg.print_ordered_edges(1,nfrac+3," 101\n");
		current_pos=bamg.print_ordered_edges(current_pos,1," 102\n");
		current_pos=bamg.print_ordered_edges(current_pos,nfrac+3," 103\n");
		bamg.print(ntotal," 1 104\n");
		bamg.close_file();

		std::ofstream fdmn( domain_filename(base_name).c_str() );
		fdmn << "EdgeDomainNames\n"
				"4\n"
				"top\n"
				"right\n"
				"bottom\n"
				"left";
		fdmn.close();
	}

};

#endif /* WAVYFRACTURE_H_ */
