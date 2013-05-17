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
	const rheolef::test v;
	const std::string domain_name;
	const double val;

public:
	NormalStressBC_RHS( const XMLConfigFile& conf, const rheolef::space& Uspace ):
		v(Uspace),
		domain_name( conf("domain_name") ),
		val( conf.atof("normal_stress_value") )
	{}

	void add_to_rhs( rheolef::field& rhs ) const
	{
		rhs += integrate( domain_name, -val*dot(v,rheolef::normal()) );
	}
};


#endif /* NORMALSTRESSBC_RHS_H_ */
