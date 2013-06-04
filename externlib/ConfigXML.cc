/*
 * ConfigXML.cc
 *
 *  Created on: 2013-04-10
 *      Author: ali
 */

#include <cstdio>
#include <stdexcept>

#include "ConfigXML.h"


XMLConfigFile::XMLConfigFile( const cstr fname ):
	root_path("/")
{
	if( !doc.LoadFile(fname) )
		throw std::runtime_error("unable to open XML file");
	rootnode = doc.RootElement();
}


XMLConfigFile::Nodeptr
XMLConfigFile::find_node( xmlpath path ) const
{
	Nodeptr nd;
	if( !node_accessible(path,nd) )
		throw std::runtime_error("XML path >>>"+root_path+path_string(path)+"<<< doesn't exist");
	return nd;
}


bool
XMLConfigFile::node_accessible( xmlpath path, Nodeptr& result ) const
{
	bool found = true;
	result = rootnode;

	for(auto name : path){
		result = result->FirstChildElement(name);
		if(!result){
			found = false;
			break;
		}
	}
	return found;
}


double
XMLConfigFile::atof_if_exist( xmlpath path, double default_val ) const
{
	return return_if_exist( path,default_val, [](cstr s){return ::atof(s);} );
}


int
XMLConfigFile::atoi_if_exist( xmlpath path, int default_val ) const
{
	return return_if_exist( path,default_val, [](cstr s){return ::atoi(s);} );
}


std::string
XMLConfigFile::path_string( xmlpath path )
{
	std::string str;
	for(const auto& s:path){
		str += "/";
		str += s;
	}
	return str;
}


//void
//XMLConfigFile::print_path( xmlpath path ) const
//{
//	printf(" >>> %s\n", (root_path+path_string(path)).c_str() );
////	for(auto name : path)
////		printf("%s/",name);
////	printf("\n");
//}
