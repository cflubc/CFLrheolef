/*
 * FlowFields.h
 *
 *  Created on: 2013-04-13
 *      Author: ali
 */

#ifndef FLOWFIELDS_H_
#define FLOWFIELDS_H_

#include <string>

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
		Qh( _omega, cf("pspace") ),
		Xh( _omega, cf("vspace"), "vector" ),
		ph(Qh, 0.),
		omega(_omega)
	{
		// Pressure is Lag multiplier we never block it's space (even if we
		// have P imposed on a boundary) or will get mass loss so blocking
		// is always for velocity space. When P is given in a boundary we add
		// pressure contribution to rhs of velocity and the correct value
		// automatically is obtained. (note this is for incompressible flow)
		bc.block_velocity_space(Xh);
		uh = rheolef::field(Xh, 0.);
		bc.set_velocity_dirichlet(uh);
	}

	void write_fields( rheolef::odiststream& o ) const
	{
		write_field(uh,"u",o);
		write_field(ph,"p",o);
	}

	geo const& get_geo() const
	{return omega;}

	std::string geo_name() const
	{return omega.name();}

	space const& Uspace() const
	{return Xh;}

	field const& Uh() const {return uh;}
	field const& Ph() const {return ph;}
	field      & Uh()       {return uh;}
	field      & Ph()       {return ph;}


	space Qh;
	space Xh;

	field ph;
	field uh;
	const geo& omega;
};



#endif /* FLOWFIELDS_H_ */
