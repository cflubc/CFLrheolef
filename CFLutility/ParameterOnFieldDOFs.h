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


/**
 * DOF naming convention for Degree Of Freedom.
 * The method begin_dof() name obtained from rheolef namin
 */
namespace ParameterOnDOFs {

using rheolef::field;
using rheolef::space;
using rheolef::Float;


template< typename T >
class unique_basic
{
	template< typename T_ >
	class unique_iterator {
	public:
		unique_iterator( T_ const& x ):
			val(x)
		{}

		void operator++() const {}

		T_ operator*() const
		{return val;}

	private:
		T_ const val;
	};

public:
	typedef unique_iterator<T> const_iterator;

	template< typename Initializer >
	unique_basic( space const& Xh, XMLConfigFile const& conf, Initializer const& I ):
		param_val( I.get(conf) )
	{}

	T get_parameter() const
	{return param_val;}

	const_iterator begin_dof() const
	{return unique_iterator<T>(param_val);}

private:
	T param_val;
};


template< typename RheoField >
class field_param_basic
{
public:
	typedef typename RheoField::const_iterator const_iterator;

	template< typename Initializer >
	field_param_basic( space const& Xh, XMLConfigFile const& conf, Initializer const& I ):
		field_param( make(Xh,conf,I) )
	{}

	RheoField get_parameter() const
	{return field_param;}

	const_iterator begin_dof() const
	{return field_param.begin_dof();}

private:
	template< typename Initializer >
	field make( space const& Xh, XMLConfigFile const& conf, Initializer const& I ){
		field f( Xh, std::numeric_limits<Float>::quiet_NaN() );
		I.get(f,conf);
		return f;
	}
	RheoField field_param;
};

typedef field_param_basic<const field> ConstNonUnique;
typedef unique_basic<const Float> ConstUnique;


////////////////////// Initializers for parameters ////////////////////////
struct unique_default_initializer {
	static Float get( XMLConfigFile const& conf )
	{return conf({},Float());}
};


struct multiRegion_field_initializer
{
	typedef std::vector<std::string> names_t;
	typedef std::vector<Float> vals_t;

	static void get( field& f, XMLConfigFile const& conf )
	{
		names_t names;
		vals_t vals;
		get_domains_and_values(conf,names,vals);
		if( names.size()!=vals.size() )
			throw std::logic_error("Number mismatch of domains/values in multiRegion parameter");
		set_field_regions_values(f,names,vals);
	}

	static void get_domains_and_values( XMLConfigFile const& conf, names_t& names, vals_t& vals ){
		conf("domain_names",&names);
		conf("values",&vals);
	}

	static void set_field_regions_values( field& f, names_t& names, vals_t& vals ){
		auto name = names.cbegin();
		for(auto val=vals.cbegin(); val!=vals.cend(); ++val,++name)
			f[*name] = *val;
	}
};

}  // namespace ParameterOnDOFs


#endif /* PARAMETERONFIELDDOFS_H_ */
