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
#include <utility>
#include <vector>

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
	typedef std::vector<double> point_arr;
	typedef char const*const tag_str;

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
	{print_edges_header(nvertices);}

	void print_edges_header( size_t const nedges )
	{println_args(file,"\nEdges ",nedges);}

	void print_subdomain_header( size_t const ndom )
	{println_args(file,"\nSubDomain ",ndom);}

	size_t print_ordered_edges( size_t const& beg,
			                  size_t const& nedges,
			                  tag_str tag ){
		size_t const last = beg+nedges;
		for(size_t i=beg; i<last; ++i)
			print(i," ",i+1,tag);
		return last;
	}

	void print_points( point_arr const& X, point_arr const& Y, tag_str tag, bool reverse=false ){
		if(!reverse){
			for(auto x=X.cbegin(), end=X.cend(), y=Y.cbegin(); x!=end; ++x,++y)
				print(*x," ",*y,tag);
		}
		else{
			for(auto x=X.crbegin(), end=X.crend(), y=Y.crbegin(); x!=end; ++x,++y)
				print(*x," ",*y,tag);
		}
	}

	size_t print_edge( size_t const& i, tag_str tag ){
		print(i," ",i+1,tag);
		return i+1;
	}

	template< typename... Args >
	void print( Args&&... args )
	{print_args(file,std::forward<Args>(args)...);}

	void close_file()
	{file.close();}
};


void
make_geo_from_bamgcad_and_dmn_file(
		std::size_t const& npoints,
		std::string const& base_name,
		std::string const& other_commandline_args );



#endif /* BAMGCAD_H_ */
