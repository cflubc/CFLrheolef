/*
 * IncompressibleNavierStokes.h
 *
 *  Created on: Aug 20, 2013
 *      Author: ali
 */

#ifndef INCOMPRESSIBLENAVIERSTOKES_H_
#define INCOMPRESSIBLENAVIERSTOKES_H_

#include <cstddef>
#include <string>
#include <iostream>
#include <cmath>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "BlockSystem_abtb.h"
#include "ConfigXML.h"
#include "ResidualTablePrinter.h"
#include "GenericIteration.h"
#include "ErrorAnalysis.h"


class IncompressibleNavierStokes
{
	typedef rheolef::Float Float;
	typedef rheolef::field field;
	typedef rheolef::space space;
	typedef rheolef::form form;
	typedef rheolef::trial trial;
	typedef rheolef::test test;
	typedef rheolef::characteristic characteristic;
	typedef rheolef::odiststream odiststream;
	typedef rheolef::quadrature_option_type quadrature;
	typedef std::string string;
	typedef std::size_t size_t;

public:
	template< typename FieldsPool, typename DirichletBC >
	IncompressibleNavierStokes( XMLConfigFile const& conf, FieldsPool& fields, DirichletBC& BC, Float const viscosity=1. ):
		Re( conf({"PhysicalParameters","Re"},Re) ),
		delta_t( conf({"PhysicalParameters","dt"},delta_t) ),
		uh(fields.Uh()),
		ph(fields.Ph()),
		dirichlet_rhs(fields.Uspace(),0.),
		uh1(uh),
		uh2(uh),
		uh_star(uh),
		v(fields.Uspace()),
		qopt( quadrature::gauss_lobatto, uh.get_space().degree() ),
		solver(conf,make_a(viscosity),make_b()),
		deltaU(uh1,uh2),
		output(30,std::cout,10,16),
		converge_loop(conf)
	{}

	void run()
	{
		solver.set_discrete_dirichlet_rhs(dirichlet_rhs,uh);
		uh1 = uh;
		converge_loop(*this);
		odiststream o (uh.get_geo().name(), "field");
		write_to_diststream(o,"u",uh, "p",ph);
		o.close();
	}

	void iterate()
	{
		uh2 = uh1;
		uh1 = uh;
		uh_star = 2.0*uh1 - uh2;
		characteristic const X1(    -delta_t*uh_star);
		characteristic const X2(-2.0*delta_t*uh_star);
//		field const l1h = integrate(dot(compose(uh1,X1),v), qopt);
//		field const l2h = integrate(dot(compose(uh2,X2),v), qopt);
		field const convec = integrate( 2.*dot(compose(uh1,X1),v) - .5*dot(compose(uh2,X2),v), qopt);
		field const lh  = dirichlet_rhs + (Re/delta_t)*convec; //(2*l1h - 0.5*l2h); // + AL.augmented_lagraniang_rhs();
		solver.solve(uh,ph,lh);
	}

	void iterate_report( size_t const niter, Float& res )
	{
		deltaU.save_field();
		iterate();
		res = deltaU.calculate_field_change()/delta_t;
//		AL.update_lagrangeMultipliers_report_strain_residual();
		output.print_header_if_needed("iteration","|Un+1-Un/dt|L2");
		output.print(niter,res);

	}

	field adapt_criteria() const
	{
		space T0h( uh.get_geo(), derivative_approx(uh.get_approx()) );
		return interpolate( T0h, sqrt(Re*norm2(uh) + 4.*norm2(D(uh))) );
	}

	form make_a( Float const& coef ) const {
		trial u( uh.get_space() );
		return integrate (2.*coef*ddot(D(u),D(v)) + 1.5*(Re/delta_t)*dot(u,v),qopt);
	}

	form make_b() const {
		trial u( uh.get_space() );
		test  q( ph.get_space() );
		return integrate (-div(u)*q, qopt);
	}

	Float const Re;
	Float const delta_t;
	field& uh;
	field& ph;
	field dirichlet_rhs;
	field uh1, uh2, uh_star;
	test const v;
	quadrature const qopt;
	BlockSystem_abtb solver;
	L2norm_calculator deltaU;
	ResidualTablePrinter<2,std::ostream> output;
	GenericIteration converge_loop;
};


#endif /* INCOMPRESSIBLENAVIERSTOKES_H_ */
