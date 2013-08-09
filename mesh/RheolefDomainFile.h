/*
 * RheolefDomainFile.h
 *
 *  Created on: Jul 21, 2013
 *      Author: ali
 */

#ifndef RHEOLEFDOMAINFILE_H_
#define RHEOLEFDOMAINFILE_H_

#include <string>
#include <fstream>
#include <initializer_list>

#include "PrintArguments.h"

class RheolefDomainFile
{
	std::ofstream file;
	typedef std::initializer_list<const char*> const name_list;

	void print_domain( std::string const& header, name_list names )
	{
		println_args(file, header,'\n', names.size() );
		for(auto& name : names)
			println_args(file,name);
		file << '\n';
	}

public:

	RheolefDomainFile( std::string const& basename ):
		file( (basename+".dmn").c_str() )
	{}

	void print_edge_domains( name_list names )
	{print_domain("EdgeDomainNames",names);}

	void print_region_domains( name_list names )
	{print_domain("RegionDomainNames",names);}

	void close_file()
	{file.close();}
};


#endif /* RHEOLEFDOMAINFILE_H_ */
