/*
 * NormalStressBC_RHS.h
 *
 *  Created on: 2013-05-17
 *      Author: ali
 */

#ifndef NORMALSTRESSBC_RHS_H_
#define NORMALSTRESSBC_RHS_H_

#include <string>
#include "rheolef.h"
#include "ConfigXML.h"

class NormalStressBC_RHS
{
	typedef rheolef::field field;

	const rheolef::test v;
	const double val;
	const field normalized_rhs;

public:
	NormalStressBC_RHS( const XMLConfigFile& conf, const rheolef::space& Uspace ):
		v(Uspace),
		val( conf.atof("normal_stress_value") ),
		normalized_rhs( -integrate(conf("edge_name"), dot(v,rheolef::normal())) )
	{}

	void add_to_rhs( field& rhs ) const
	{rhs += get_rhs();}

	field get_rhs( double const& x ) const
	{return x*normalized_rhs;}

	field get_rhs() const
	{return val*normalized_rhs;}
};


#endif /* NORMALSTRESSBC_RHS_H_ */
