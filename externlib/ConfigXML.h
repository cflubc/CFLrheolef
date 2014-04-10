/**
 * @file SimpleXML.h
 *
 * @date Jun 14, 2011
 * @author ali
 */

#ifndef CONFIGXML_H_
#define CONFIGXML_H_

#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <initializer_list>

#include "tinyxml.h"


class XMLConfigFile
{
	typedef char const* cstr;
    typedef std::initializer_list<cstr> xmlpath;
    typedef TiXmlElement const* Nodeptr;

    std::string const root_path;
    Nodeptr rootnode;
    TiXmlDocument doc;

    cstr get_txt( xmlpath path ) const;
    static std::string path_string( xmlpath path );

public:
    bool node_accessible( xmlpath path, Nodeptr& result ) const;
    Nodeptr find_node( xmlpath path ) const;

    XMLConfigFile( Nodeptr node, std::string const& parent_path ):
    	root_path(parent_path),
    	rootnode(node)
    {}

    XMLConfigFile( const cstr fname );
    XMLConfigFile child( xmlpath path ) const;

    cstr operator()( xmlpath  path ) const;
    double atof( xmlpath path ) const;
    int    atoi( xmlpath path ) const;

    double atof_if_exist( xmlpath path, double default_val ) const;
    int    atoi_if_exist( xmlpath path, int    default_val ) const;

    template< typename T, typename Function >
    T return_if_exist( xmlpath path, T default_val, const Function& f ) const;

    cstr return_txt_if_exist( xmlpath path, cstr const default_val ) const;

	template< typename T >
	void operator()( xmlpath path, T *const val ) const;

	template< typename T >
	T operator()( xmlpath path, T val ) const;

	/**
	 * This overload is necessary: e.g. </opt>-hmax .5</opt> and we want to
	 * get one string "-hmax .5". If use the normal templated function then
	 * we only read the first "-hmax" from the stringstream.
	 */
	void operator()( xmlpath path, std::string *const val ) const;

    template< typename T >
    void set_if_path_exist( xmlpath path, T *const val ) const;

    template< typename T >
    T get_if_path_exist( xmlpath path, T default_val ) const;


    ///< convinient functions to avoid writing {} when path is only one element
    XMLConfigFile child( cstr one_path ) const;
    cstr operator()( cstr one_path ) const;
    double atof( cstr one_path ) const;
    int    atoi( cstr one_path ) const;
    double atof_if_exist( cstr one_path, double default_val ) const;
    int    atoi_if_exist( cstr one_path, int    default_val ) const;

    template< typename T >
	void operator()( cstr one_path, T *const val ) const;

    template< typename T>
    T operator()( cstr one_path, T val ) const;

    template< typename T >
    void set_if_path_exist( cstr one_path, T *const val ) const;

    template< typename T >
    T get_if_path_exist( cstr one_path, T default_val ) const;

};



inline
XMLConfigFile::cstr
XMLConfigFile::get_txt( xmlpath path ) const
{return find_node(path)->GetText();}

inline
XMLConfigFile::cstr
XMLConfigFile::operator()( xmlpath path ) const
{return get_txt(path);}

inline
XMLConfigFile
XMLConfigFile::child( xmlpath path ) const
{ return XMLConfigFile( find_node(path), root_path+path_string(path) ); }


template< typename T >
void XMLConfigFile::operator()( xmlpath path, T *const val ) const
{
	std::istringstream ss( get_txt(path) );
	ss >> *val;
}


inline
void XMLConfigFile::operator()( xmlpath path, std::string *const val ) const
{*val = operator()(path);}


template< typename T >
T XMLConfigFile::operator()( xmlpath path, T ) const
{
	T val;
	operator()(path,&val);
	return val;
}


template< typename T >
void XMLConfigFile::set_if_path_exist( xmlpath path, T *const val ) const
{
	Nodeptr nd;
	if( node_accessible(path,nd) )
		operator()(path,val);
}


template< typename T >
T XMLConfigFile::get_if_path_exist( xmlpath path, T default_val ) const
{
	Nodeptr nd;
	if( node_accessible(path,nd) )
		operator()(path,&default_val);
	return default_val;
}


inline int
XMLConfigFile::atoi( xmlpath path ) const
{return ::atoi( get_txt(path) );}

inline double
XMLConfigFile::atof( xmlpath path ) const
{return ::atof( get_txt(path) );}


template< typename T, typename Function >
T XMLConfigFile::return_if_exist( xmlpath path, T default_val, const Function& f ) const
{
	Nodeptr nd;
	if( node_accessible(path,nd) )
		return f( nd->GetText() );
	return default_val;
}

inline
XMLConfigFile::cstr
XMLConfigFile::return_txt_if_exist( xmlpath path, cstr const default_val ) const
{
	// using c++11 lambda function
	return return_if_exist(path,default_val,[](cstr s){return s;});
}


template< typename T >
void operator>>( std::istream& is, std::vector<T>& vec )
{
	while( !is.eof() )
	{
		T tmp;
		is >> tmp;
		vec.push_back(tmp);
	}
}

// to convert cstr to initializer_list with one element it is necessary to
// write type explicity as xmlpath{path}, otherwise {path} is still
// interpreted as cstr and infinite recursive template instantiation happens!
//
// auto x = {2};   // x: initializer_list of int, OK
// int  x = {2};   // x: an int!, OK

inline XMLConfigFile::cstr
XMLConfigFile::operator()( cstr one_path ) const
{return get_txt( xmlpath{one_path} );}

inline XMLConfigFile
XMLConfigFile::child( cstr one_path ) const
{return child( xmlpath{one_path} );}

template< typename T >
inline void
XMLConfigFile::operator()( cstr one_path, T *const val ) const
{operator()(xmlpath{one_path},val);}

template< typename T>
inline T
XMLConfigFile::operator()( cstr one_path, T val ) const
{return operator()(xmlpath{one_path},val);}


template< typename T >
inline T
XMLConfigFile::get_if_path_exist( cstr one_path, T default_val ) const
{return get_if_path_exist(xmlpath{one_path},default_val);}


inline
int XMLConfigFile::atoi( cstr one_path ) const
{ return atoi(xmlpath{one_path}); }

inline
double XMLConfigFile::atof( cstr one_path ) const
{ return atof(xmlpath{one_path}); }

inline
double XMLConfigFile::atof_if_exist( cstr one_path, double default_val ) const
{ return atof_if_exist(xmlpath{one_path},default_val); }

inline
int XMLConfigFile::atoi_if_exist( cstr one_path, int default_val ) const
{ return atoi_if_exist(xmlpath{one_path},default_val); }


#define XML_VAL(xmlconf,var,path) \
	xmlconf(path,var)

#define XML_VAL_IF_EXIST(xmlconf,default_val,path) \
	xmlconf.get_if_path_exist({path},default_val)


#define XML_INIT_VAR(xmlconf,var,path) \
	var( XML_VAL(xmlconf,var,path) )

#define XML_INIT_VAR_IF_EXIST(xmlconf,var,default_val,path) \
	var( XML_VAL_IF_EXIST(xmlconf,default_val,path) )


#endif

