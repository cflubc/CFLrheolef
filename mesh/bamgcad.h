/*
 * bamgcad.h
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */

#ifndef BAMGCAD_H_
#define BAMGCAD_H_

#include <cstdlib>
#include <fstream>
#include <string>

#include "PrintArguments.h"


inline std::string
bamgcad_filename( std::string const& base )
{return base+".bamgcad";}

inline std::string
bamgmesh_filename( std::string const base )
{return base+".bamg";}

/**
 * Helper class for generating bamg geometry input file
 */
class bamgcad
{
	typedef std::string string;

	size_t const nvertices;
	std::ofstream file;

public:

	bamgcad( size_t const n, std::string const& base_name ):
		nvertices(n),
		file( bamgcad_filename(base_name).c_str() )
	{
		println_args(file,
				"MeshVersionFormatted 0\n"
				"AngleOfCornerBound 0\n"
				"Dimension 2\n\n"
				"Vertices ", nvertices);
	}

	void print_edges_header()
	{println_args(file,"\nEdges ",nvertices);}

	void print_ordered_edges( size_t const& beg,
			                  size_t const& nedges,
			                  char const*const tag ){
		for(size_t i=beg; i<beg+nedges; ++i)
			print_args(file,i," ",i+1,tag);
	}

	template< typename... Args >
	void print( Args const&... args )
	{print_args(file,args...);}

	void close_file()
	{file.close();}
};


void
make_geo_from_bamgcad_and_dmn_file(
		std::string const& base_name,
		std::string const& other_commandline_args );



#endif /* BAMGCAD_H_ */
