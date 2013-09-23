/*
 * ErrorAnalysis.h
 *
 *  Created on: 2013-04-22
 *      Author: ali
 */

#ifndef ERRORANALYSIS_H_
#define ERRORANALYSIS_H_

#include "rheolef.h"
#include "CFL.h"


class L2norm_calculator
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

public:

	L2norm_calculator( field const& _f );
	L2norm_calculator( field const& f1, field const& f2 );

	void save_field();

	/// This function is invoked after f has changed
	rheolef::Float calculate_field_change();

	rheolef::Float calculate_fieldL2( field const& q ) const
	{return rheolef::sqrt(mass(q,q));}

private:

	bool const initialized_with_one_field;
	field const& f;
	field const*const f_new;
	field field_delta;
	rheolef::form const mass;
};


#endif /* ERRORANALYSIS_H_ */
