/*
 * FlowFields.h
 *
 *  Created on: 2013-04-13
 *      Author: ali
 */

#ifndef FLOWFIELDS_H_
#define FLOWFIELDS_H_


#include "rheolef/diststream.h"
#include "rheolef.h"

#include "CFL.h"
#include "ConfigXML.h"


struct FlowFields
{
	typedef rheolef::space space;
	typedef rheolef::field field;
	typedef rheolef::geo     geo;

	template< typename VelDirichletBC >
	FlowFields( const XMLConfigFile& cf , const geo& _omega, const VelDirichletBC& bc ):
		Qh( _omega, cf({"FEfields","pspace"}) ),
		Xh( _omega, cf({"FEfields","vspace"}), "vector" ),
		Ph(Qh, 0.),
		omega(_omega)
	{
		// Pressure is Lag multiplier we never block it's space (even if we
		// have P imposed on a boundary) or will get mass loss so blocking
		// is always for velocity space. When P is given in a boundary we add
		// pressure contribution to rhs of velocity and the correct value
		// automatically is obtained. (note this is for incompressible flow)
		bc.block_velocity_space(Xh);
		Uh = rheolef::field(Xh, 0.);
		bc.set_velocity_dirichlet(Uh);
	}

	void write_fields( rheolef::odiststream& o ) const
	{
		write_field(Uh,"u",o);
		write_field(Ph,"p",o);
	}

	space Qh;
	space Xh;

	field Ph;
	field Uh;
	const geo& omega;
};



#endif /* FLOWFIELDS_H_ */
