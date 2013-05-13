/*
 * ErrorAnalysis.cc
 *
 *  Created on: 2013-05-05
 *      Author: ali
 */

#include "ErrorAnalysis.h"


rheolef::Float
L2norm_calculator::calculate_field_change()
{
	field_old -= *f;
	rheolef::Float norm = rheolef::sqrt( mass(field_old,field_old) );
	field_old = *f;
	return norm;
}



