/*
 * AdaptationCriterions.cc
 *
 *  Created on: 2013-04-16
 *      Author: ali
 */

#include "CFL.h"
#include "adaptationCriterions.h"


rheolef::field stokes_criterion(const rheolef::field& uh)
{
	using namespace rheolef;
	const std::string approx = derivative_approx( uh.get_approx() );
	space Th( uh.get_geo(), approx, "tensor" );
	form _2D( uh.get_space(), Th, "2D" );
	form inv_mt( Th, Th, "inv_mass" );
	field _2Duh  = inv_mt*(_2D*uh);

	space X0h( uh.get_geo(), approx );
	return interpolate( X0h, sqrt(norm2(_2Duh)) );
}
