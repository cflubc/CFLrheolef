/*
 * ErrorAnalysis.cc
 *
 *  Created on: 2013-05-05
 *      Author: ali
 */

#include "ErrorAnalysis.h"


L2norm_calculator::L2norm_calculator( field const& _f ):
	initialized_with_one_field(true),
	f(_f),
	f_new(nullptr),
	field_delta(_f),
	mass(_f.get_space(),_f.get_space(),"mass")
{}


L2norm_calculator::L2norm_calculator( field const& f1, field const& f2 ):
	initialized_with_one_field(false),
	f(f1),
	f_new(&f2),
	mass(f1.get_space(),f1.get_space(),"mass")
{}


void
L2norm_calculator::save_field()
{
	if( initialized_with_one_field )
		field_delta = f;
}


rheolef::Float
L2norm_calculator::calculate_field_change()
{
	if( initialized_with_one_field )
		field_delta -= f;
	else
		field_delta = *f_new - f;
	return calculate_fieldL2(field_delta);
}



