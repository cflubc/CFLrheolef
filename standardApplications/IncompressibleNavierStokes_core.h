/*
 * IncompressibleNavierStokes_core.h
 *
 *  Created on: Oct 8, 2013
 *      Author: ali
 */

#ifndef INCOMPRESSIBLENAVIERSTOKES_CORE_H_
#define INCOMPRESSIBLENAVIERSTOKES_CORE_H_


#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "BlockSystem_abtb.h"


class IncompressibleNavierStokes_core
{
public:

	typedef rheolef::Float Float;
	typedef rheolef::field field;
	typedef rheolef::space space;
	typedef rheolef::form form;
	typedef rheolef::trial trial;
	typedef rheolef::test test;
	typedef rheolef::characteristic characteristic;
	typedef rheolef::quadrature_option_type quadrature;
	typedef std::string string;
	typedef std::size_t size_t;

	template< typename FieldsPool, typename DirichletBC >
	IncompressibleNavierStokes_core( XMLConfigFile const& conf, FieldsPool& fields, DirichletBC& BC, Float const viscosity ):
		Re( conf({"PhysicalParameters","Re"},Re) ),
		delta_t( conf({"PhysicalParameters","dt"},delta_t) ),
		uh(fields.Uh()),
		ph(fields.Ph()),
		dirichlet_rhs(fields.Uspace(),0.),
		convect_rhs(dirichlet_rhs),
		uh1(uh),
		uh2(uh),
		uh_star(uh),
		v(fields.Uspace()),
		qopt( quadrature::gauss_lobatto, uh.get_space().degree() ),
		solver(conf,make_a(viscosity),make_b())
	{}


	void solve( field const& rhs )
	{solver.solve(uh,ph,rhs);}

	void compute_convective_rhs()
	{
		uh2 = uh1;
		uh1 = uh;
		uh_star = 2.0*uh1 - uh2;
		characteristic const X1(    -delta_t*uh_star);
		characteristic const X2(-2.0*delta_t*uh_star);
		convect_rhs = (Re/delta_t)*integrate( 2.*dot(compose(uh1,X1),v) - .5*dot(compose(uh2,X2),v), qopt);
	}

	void write( rheolef::odiststream& o ) const
	{write_to_diststream(o,"u",uh, "p",ph);}

	Float const Re;
	Float const delta_t;
	field& uh;
	field& ph;
	field dirichlet_rhs;
	field convect_rhs;
	field uh1, uh2, uh_star;
	test const v;
	quadrature const qopt;
	BlockSystem_abtb solver;

private:

	form make_a( Float const& coef ) const {
		trial u( uh.get_space() );
		return integrate (2.*coef*ddot(D(u),D(v)) + 1.5*(Re/delta_t)*dot(u,v),qopt);
	}

	form make_b() const {
		trial u( uh.get_space() );
		test  q( ph.get_space() );
		return integrate (-div(u)*q, qopt);
	}
};



#endif /* INCOMPRESSIBLENAVIERSTOKES_CORE_H_ */
