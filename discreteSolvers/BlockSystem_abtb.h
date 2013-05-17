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
#include <cassert>

typedef rheolef::vec<rheolef::Float> rheo_vec;

class abtb_solver_preconditionedP
{
	typedef rheolef::form form;

	form mp;
	rheolef::solver_abtb solver;

public:
	abtb_solver_preconditionedP( const XMLConfigFile& cf, const form& a, const form& b ):
		mp(b.get_second_space(),b.get_second_space(),"mass"),
		solver( a.uu(), b.uu(), mp.uu() )
	{}

	void solve( const rheo_vec& urhs, const rheo_vec& prhs, rheo_vec& u, rheo_vec& p ) const
	{solver.solve(urhs,prhs,u,p);}
};


class BlockSystem_abtb
{
	typedef rheolef::form form;
	typedef rheolef::field field;

	form a;
	form b;
	rheo_vec p_rhs;
	abtb_solver_preconditionedP abtb_method;

public:
	BlockSystem_abtb( const XMLConfigFile& cf, const form& _a, const form& _b ):
		a(_a),
		b(_b),
		abtb_method(cf,a,b)
	{}

	void solve( field& Uh, field& Ph, field& rhs ) const
	{
		abtb_method.solve( rhs.u(), p_rhs, Uh.set_u(), Ph.set_u() );
	}

	void set_discrete_dirichlet_rhs( field& urhs, const field& Uh ){
		urhs.set_u() = -a.ub()*Uh.b();
		p_rhs = -( b.ub()*Uh.b() );
	}
};


#endif /* PCG_ABTB_H_ */
