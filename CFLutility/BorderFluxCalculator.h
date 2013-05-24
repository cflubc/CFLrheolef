/*
 * BorderFluxCalculator.h
 *
 *  Created on: 2013-05-21
 *      Author: ali
 */

#ifndef BORDERFLUXCALCULATOR_H_
#define BORDERFLUXCALCULATOR_H_

#include <string>
#include "rheolef.h"
#include "CFL.h"

class BorderFluxCalculator
{
	typedef rheolef::field field;

	std::string const border_name;
	field const& f;
	field const m;

	field vdotn( std::string const& name ) const {
		rheolef::test v( f.get_space() );
		return integrate( border_name, dot(rheolef::normal(),v) );
	}

public:
	BorderFluxCalculator( std::string const& name, field const& _f ):
		border_name(name),
		f(_f),
		m( vdotn(border_name) )
	{}

	rheolef::Float calc_flux() const {
		return vector_dot( m[border_name], f[border_name] );
	}
};



#endif /* BORDERFLUXCALCULATOR_H_ */
