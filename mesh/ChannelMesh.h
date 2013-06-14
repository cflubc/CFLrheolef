/*
 * ChannelMesh.h
 *
 *  Created on: 2013-05-27
 *      Author: ali
 */

#ifndef CHANNELMESH_H_
#define CHANNELMESH_H_


#include <string>
#include <fstream>

#include "ConfigXML.h"
#include "bamgcad.h"


class ChannelMesh
{
public:

	ChannelMesh( XMLConfigFile const& conf, std::string const& base_name ):
		whalf( .5*conf.atof("width") ),
		lhalf( .5*conf.atof("length") )
	{
		bamgcad bamg(4,base_name);
		std::string const type( conf("type") );
		if( "full"==type )
			gen_rectangle( -lhalf,-whalf, lhalf,whalf, bamg );
		else if( "symy"==type )
			gen_rectangle( -lhalf,0, lhalf,whalf, bamg );
		else if( "symx"==type )
			gen_rectangle( 0,-whalf, lhalf,whalf, bamg );
		else if( "symxy"==type )
			gen_rectangle( 0,0, lhalf,whalf, bamg );
		else
			throw std::logic_error("Wrong type for ChannelMesh provided");

		std::ofstream fdmn( domain_filename(base_name) );
		fdmn <<  "EdgeDomainNames\n"
		         "4\n"
		         "bottom\n"
		         "right\n"
		         "top\n"
		         "left";
		fdmn.close();
	}

private:

	double const whalf;
	double const lhalf;

	static void gen_rectangle( double const& x_bl, double const& y_bl,
						double const& x_tr, double const& y_tr, bamgcad& bamg )
	{
		bamg.print(x_bl," ",y_bl," 1\n");
		bamg.print(x_tr," ",y_bl," 2\n");
		bamg.print(x_tr," ",y_tr," 3\n");
		bamg.print(x_bl," ",y_tr," 4\n");

		bamg.print_edges_header();
		bamg.print(
				"1 2  101\n"
				"2 3  102\n"
				"3 4  103\n"
				"4 1  104\n" );
		bamg.close_file();
	}
};


#endif /* CHANNELMESH_H_ */
