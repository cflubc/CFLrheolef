/*
 * ParameterOnFieldDOFs.h
 *
 *  Created on: Jul 20, 2013
 *      Author: ali
 */

#ifndef PARAMETERONFIELDDOFS_H_
#define PARAMETERONFIELDDOFS_H_


#include <limits>
#include <vector>
#include <string>
#include <stdexcept>

#include "rheolef.h"
#include "ConfigXML.h"


/** DOF naming convention for Degree Of Freedom */
namespace ParameterOnFieldDOFs {

class Unique
{
	typedef rheolef::Float Float;
	Float const p;

public:
	Unique( rheolef::space const& Xh, XMLConfigFile const& conf ):
		p( conf({},p) )
	{}

	void begin() const {}
	void operator++() const {}
	Float val_on_current_dof() const {return p;}
	Float parameter() const {return p;}
};


template< typename initializer >
class Constant
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

	field const f;
	field::const_iterator iter;

	field init(rheolef::space const& Xh, XMLConfigFile const& conf ){
		field f(Xh, std::numeric_limits<Float>::quiet_NaN());
		initializer::init_field(f,conf);
		return f;
	}

public:
	Constant( rheolef::space const& Xh, XMLConfigFile const& conf ):
		f( init(Xh,conf) )
	{}

	void begin() {iter = f.begin_dof();}
	void operator++()  {++iter;}
	Float val_on_current_dof() const  {return *iter;}
	field parameter() const {return f;}
};


struct multi_region_initializer
{
	static void init_field( rheolef::field& f, XMLConfigFile const& conf )
	{
		using namespace std;
		vector<string> names;
		vector<rheolef::Float> vals;
		conf("domain_names",&names);
		conf("values",&vals);
		if( names.size()!=vals.size() )
			throw logic_error("Number of domains is different from values");

		auto name = names.cbegin();
		for(auto val=vals.cbegin(); val!=vals.cend(); ++val,++name)
			f[*name] = *val;
	}
};

typedef Constant<multi_region_initializer> MultiRegion;
}  // namespace ParameterOnFieldDOFs


#endif /* PARAMETERONFIELDDOFS_H_ */
