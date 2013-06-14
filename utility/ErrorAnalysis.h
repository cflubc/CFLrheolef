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


struct L2norm_calculator
{
	typedef rheolef::field field;

	field const* f;
	field field_old;
	rheolef::form mass;

	L2norm_calculator( const field& _f ):
		f(&_f),
		mass(_f.get_space(),_f.get_space(),"mass")
	{}

	void save_field()
	{field_old = *f;}

	rheolef::Float calculate_fieldL2( field const& f ) const
	{return rheolef::sqrt(mass(f,f));}

	/// This function is invoked after f has changed
	rheolef::Float calculate_field_change();
};



#endif /* ERRORANALYSIS_H_ */
