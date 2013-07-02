/*
 * BodyForce.h
 *
 *  Created on: 2013-06-02
 *      Author: ali
 */

#ifndef BODYFORCE_H_
#define BODYFORCE_H_

#include "rheolef.h"

#include "CFL.h"
#include "ConfigXML.h"


class BodyForce
{
	typedef rheolef::field field;

	const rheolef::test v;
	const rheolef::point bodyf_vector;
	const field normalized_rhs;

public:
	enum : bool { isLinear=true };

	BodyForce( XMLConfigFile const& conf, rheolef::space const& Uspace ):
		v(Uspace),
		XML_INIT_VAR(conf,bodyf_vector,"bodyforce_vector"),
		normalized_rhs( integrate(rheolef::dot(bodyf_vector,v)) )
	{}

	void add_to_rhs( field& rhs ) const
	{rhs += get_rhs();}

	field get_rhs( double const& scale ) const
	{return scale*normalized_rhs;}

	field get_rhs() const
	{return normalized_rhs;}
};


#endif /* BODYFORCE_H_ */
