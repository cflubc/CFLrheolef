/*
 * BCs.h
 *
 *  Created on: 2013-04-13
 *      Author: ali
 */

#ifndef BCS_H_
#define BCS_H_

#include "rheolef.h"
#include "ConfigXML.h"

struct cavityBC
{
	void block_velocity_space( rheolef::space& Xh ) const {
		Xh.block("top");
		Xh.block("bottom");
		Xh.block("left");
		Xh.block("right");
	}

	void set_velocity_dirichlet( rheolef::field& Uh ) const
	{Uh[0]["top"] = 1.;}

	cavityBC( XMLConfigFile const& conf )
	{}
};


struct channel_fullBC
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

	channel_fullBC( XMLConfigFile const& conf )
	{}
};


struct bubble_BC
{
	void block_velocity_space( rheolef::space& Xh ) const {
		Xh.block("top");
		Xh.block("bottom");
		Xh[1].block("left");
		Xh[1].block("right_top");
		Xh[1].block("right_bottom");
	}

	void set_velocity_dirichlet( rheolef::field& Uh ) const {
		Uh["top"] = 0.;
		Uh["bottom"] = 0.;
		Uh[1]["left"] = 0.;
		Uh[1]["right_top"] = 0.;
		Uh[1]["right_bottom"] = 0.;
	}

	bubble_BC( XMLConfigFile const& )
	{}
};


#endif /* BCS_H_ */


