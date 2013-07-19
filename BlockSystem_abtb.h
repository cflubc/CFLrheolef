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

#include <cassert>
#include <cstddef>
#include <iostream>

#include "rheolef.h"
#include "rheolef/solver_abtb.h"

#include "ConfigXML.h"


typedef rheolef::vec<rheolef::Float> rheo_vec;

template<
    class MatrixSolver,
    class Matrix,
    class Vector,
    class Real,
    class Size>
int
uz_abtb(
    const MatrixSolver&          m_solver,
    const Matrix&                b,
    Vector&                      u,
    Vector&                      p,
    const Vector&                f,
    const Vector&                g,
    const Real&                  rho,
    Size&                        max_iter,
    Real&                        tol,
    std::ostream*                p_cres = 0)
{
    Real residu = 1.;
    for (Size k=1; k <= max_iter; k++) {
        u = m_solver.solve(f - b.trans_mult(p));
        Vector R = b*u - g;
        p += rho*R ;
        residu = norm(R);
        if (p_cres) *p_cres << "[uzawa_abtb] " << k << " " << residu << "\n" ;
        if (residu <= tol) {
			u = m_solver.solve(f - b.trans_mult(p));
			tol = residu;
			max_iter = k;
			return 0;
        }
    }
    tol = residu;
    return 1;
}

class uzawa_abtb_solver
{
	typedef rheolef::Float Float;
	typedef rheolef::form form;

	Float const rho;
	Float const tolerance;
	std::size_t const max_iter;
	form const& b;
	form const ar;
	rheolef::solver const svr;

	form augmented_form( form const& a, form const& b ) const{
		form ar = a+rho*trans(b)*b;
		ar.uu().set_symmetry(true);
		return ar;
	}

public:

	uzawa_abtb_solver( XMLConfigFile const& cf, form const& a, form const& _b ):
//		XML_INIT_VAR(cf,rho,"rho"),
//		XML_INIT_VAR(cf,tolerance,"tolerance"),
//		XML_INIT_VAR(cf,max_iter,"max_iteration"),
		rho(1e7),
		tolerance(1e-11),
		max_iter(30),
		b(_b),
		ar( augmented_form(a,_b) ),
		svr( ldlt(ar.uu()) )
	{}

	void solve( const rheo_vec& urhs, const rheo_vec& prhs, rheo_vec& u, rheo_vec& p ) const
	{
		Float tol = tolerance;
		Float maxitr = max_iter;
		uz_abtb(svr,b.uu(),u,p,urhs,prhs,rho,maxitr,tol);
	}
};






class abtb_solver_preconditionedP
{
	typedef rheolef::form form;

	form const mp;
	rheolef::solver_abtb const solver;

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

	void solve( field& Uh, field& Ph, field const& rhs ) const
	{
		abtb_method.solve( rhs.u(), p_rhs, Uh.set_u(), Ph.set_u() );
	}

	void set_discrete_dirichlet_rhs( field& urhs, const field& Uh ){
		urhs.set_b() = 0.;
		urhs.set_u() = -a.ub()*Uh.b();
		p_rhs = -( b.ub()*Uh.b() );
	}
};


#endif /* PCG_ABTB_H_ */
