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

	cavityBC( const XMLConfigFile& conf )
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
		Uh[0]["left"] = 0.;
		Uh[0]["right"] = 0.;
	}

	channel_fullBC( const XMLConfigFile& conf )
	{}
};

#endif /* BCS_H_ */
