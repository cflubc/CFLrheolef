/**
 * @file pCG_abtb.h
 *  A wrapper class providing a unified interface for
 *  solver_abtb of rheolef
 *
 * @date 2013-04-13
 * @author ali
 */

#ifndef PCG_ABTB_H_
#define PCG_ABTB_H_

#include "rheolef.h"
#include "rheolef/solver_abtb.h"
#include "ConfigXML.h"


class BlockSystem_abtb
{
	typedef rheolef::form form;
	typedef rheolef::field field;

	form a;
	form b;
	form mp;
	rheolef::solver_abtb sol;

	void slv( field& Uh, field& Ph, decltype(Uh.u()) rhs ) const
	{sol.solve( rhs, -(b.ub()*Uh.b()), Uh.set_u(), Ph.set_u() );}

public:
	BlockSystem_abtb( const XMLConfigFile& cf, const form& _a, const form& _b ):
		a(_a),
		b(_b),
		mp(_b.get_second_space(),_b.get_second_space(),"mass"),
		sol( a.uu(), b.uu(), mp.uu() )
	{}

	void solve( field& Uh, field& Ph ) const
	{slv( Uh, Ph, -(a.ub()*Uh.b()) );}

	void solve( field& Uh, field& Ph, const field& rhs ) const
	{slv( Uh, Ph, rhs.u()-(a.ub()*Uh.b()) );}
};


#endif /* PCG_ABTB_H_ */
