/*
 * BCs.h
 *
 *  Created on: 2013-04-13
 *      Author: ali
 */

#ifndef BCS_H_
#define BCS_H_

#include <stdexcept>
#include "rheolef.h"

#include "ConfigXML.h"
#include "CFL.h"


struct cavityBC
{
	void block_velocity_space( rheolef::space& Xh ) const {
		Xh.block("top");
		Xh.block("bottom");
		Xh.block("left");
		Xh.block("right");
	}

	void set_velocity_dirichlet( rheolef::field& Uh ) const
	{
		Uh["left"] = 0.;
		Uh["right"] = 0.;
		Uh[0]["top"] = 1.;
		Uh[1]["top"] = 0.;
		Uh["bottom"] = 0.;
	}

	cavityBC( XMLConfigFile const& conf )
	{}
};


class channelBC
{
	typedef rheolef::space space;
	typedef rheolef::field field;

public:

	void block_velocity_space( space& Xh ) const {
		Xh.block("top");
		Xh[1].block("left");
		Xh[1].block("right");
		switch (type) {
			case sym:
				Xh[1].block("bottom");
				break;

			case full:
				Xh.block("bottom");
				break;
		}

	}

	void set_velocity_dirichlet( field& Uh ) const {
		Uh["top"] = 0.;
		Uh[1]["left"] = 0.;
		Uh[1]["right"] = 0.;
		switch (type) {
			case sym:
				Uh[1]["bottom"] = 0.;
				break;

			case full:
				Uh["bottom"] = 0.;
				break;
		}
	}

	channelBC( XMLConfigFile const& conf )
	{
		cstr const t = conf("BCtype");
		if( cstrcmp(t,"symmetry") )
			type = sym;
		else if( cstrcmp(t,"full") )
			type = full;
		else
			throw std::logic_error("Wrong type for channelBC, can only be symmetry/full");
	}

	enum BCType:int { sym, full };
	BCType type;
};


struct bubble_BC
{
	void block_velocity_space( rheolef::space& Xh ) const {
		Xh.block("top");
		Xh.block("bottom");
		Xh[1].block("left");
		Xh[1].block("right");
	}

	void set_velocity_dirichlet( rheolef::field& Uh ) const {
		Uh["top"] = 0.;
		Uh["bottom"] = 0.;
		Uh[1]["left"] = 0.;
		Uh[1]["right"] = 0.;
	}

	bubble_BC( XMLConfigFile const& )
	{}
};


#endif /* BCS_H_ */


