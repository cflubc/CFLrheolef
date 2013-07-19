/*
 * NormalStressBC_RHS.h
 *
 *  Created on: 2013-05-17
 *      Author: ali
 */

#ifndef NORMALSTRESSBC_RHS_H_
#define NORMALSTRESSBC_RHS_H_

#include "rheolef.h"
#include "ConfigXML.h"


class NormalStressBC_RHS
{
	typedef rheolef::field field;

	const rheolef::test v;
	const double val;
	const field normalized_rhs;
	double scale;

public:
	enum : bool { isLinear=true };

	NormalStressBC_RHS( XMLConfigFile const& conf, rheolef::space const& Uspace ):
		v(Uspace),
		XML_INIT_VAR_IF_EXIST(conf,val,0.,"normal_stress_value"),
		normalized_rhs( -integrate(conf("edge_name"), dot(v,rheolef::normal())) ),
		scale(1.)
	{}

	void add_to_rhs( field& rhs ) const
	{rhs += scale*get_rhs();}

	field get_rhs( double const& x ) const
	{return x*normalized_rhs;}

	field get_rhs() const
	{return val*normalized_rhs;}

	void set_scale_factor( double const& x )
	{scale=x;}
};


#endif /* NORMALSTRESSBC_RHS_H_ */
