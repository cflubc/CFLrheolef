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

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "BlockSystem_abtb.h"
#include "ConfigXML.h"
#include "ResidualTablePrinter.h"

template< typename BasicAugmentedLagrangian >
class IncompressibleNavierStokes
{
	typedef rheolef::Float Float;
	typedef rheolef::field field;
	typedef rheolef::space space;
	typedef rheolef::form form;
	typedef rheolef::characteristic characteristic;
	typedef rheolef::odiststream odiststream;
	typedef std::string string;
	typedef std::size_t size_t;

	void navier_stokes_solve()
	{
	  const space& Xh = uh.get_space();
	  const space& Qh = ph.get_space();
	  rheolef::quadrature_option_type qopt;
	  qopt.set_family(rheolef::quadrature_option_type::gauss_lobatto);
	  qopt.set_order(Xh.degree());
	  rheolef::trial u(Xh);
	  rheolef::test  v(Xh), q(Qh);
	  form const m  = integrate (dot(u,v), qopt);
	  form const a  = integrate (10.*2.*ddot(D(u),D(v)) + 1.5*(Re/delta_t)*dot(u,v), qopt);
	  form const b  = integrate (-div(u)*q, qopt);
	  BlockSystem_abtb solver(xml,a,b);
	  field dirichlet_rhs(Xh,0.);
	  solver.set_discrete_dirichlet_rhs(dirichlet_rhs,uh);

	  auto output = make_residual_table(25,std::cout,10,16);
	  bool converged = false;
	  field uh1 = uh;
	  size_t n = 0;
	  Float residual;
	  do
	  {
		field uh2;
		for(size_t i=0; i<report_freq; ++i)
		{
			uh2 = uh1;
			uh1  = uh;
			field const uh_star = 2.0*uh1 - uh2;
			characteristic const X1(    -delta_t*uh_star);
			characteristic const X2(-2.0*delta_t*uh_star);
			field const l1h = integrate(dot(compose(uh1,X1),v), qopt);
			field const l2h = integrate(dot(compose(uh2,X2),v), qopt);


			AL.update_lagrangeMultipliers_fast();
			field const lh  = dirichlet_rhs + (Re/delta_t)*(2*l1h - 0.5*l2h) + AL.augmented_lagraniang_rhs();
			solver.solve(uh,ph,lh);
		}
		n += report_freq;

	    field const duh_dt = (3.*uh - 4.*uh1 + uh2)/(2.*delta_t);
	    residual = sqrt(m(duh_dt,duh_dt));
	    output.print_header_if_needed("\niteration","|du/dt|L2");
	    output.print(n,residual);
	  } while( n<max_iter && tol<residual );
	  print_solution_convergence_message(converged);
	}
public:
	enum : bool { isLinear=false };

	template< typename FieldsPool, typename DirichletBC >
	IncompressibleNavierStokes( XMLConfigFile const& conf, FieldsPool& fields, DirichletBC& BC ):
	uh(fields.Uh()),
	ph(fields.Ph()),
	Re( conf({"PhysicalParameters","Re"},Re) ),
	delta_t( conf({"PhysicalParameters","dt"},delta_t) ),
	XML_INIT_VAR(conf,tol,"convergence_limit"),
	XML_INIT_VAR(conf,max_iter,"max_iteration"),
	XML_INIT_VAR(conf,report_freq,"report_frequency"),
	xml(conf),
	AL(conf,fields,BC)
	{}

	void run()
	{
		navier_stokes_solve();
		odiststream o (uh.get_geo().name(), "field");
		write_to_diststream(o, "Re",Re, "delta_t",delta_t, "u",uh, "p",ph);
		o.close();
	}

	field adapt_criteria() const
	{}

	field& uh;
	field& ph;
	Float const Re;
	Float const delta_t;
	Float const tol;
	size_t const max_iter;
	size_t const report_freq;
	XMLConfigFile xml;
	BasicAugmentedLagrangian AL;
};


#endif /* INCOMPRESSIBLENAVIERSTOKES_H_ */
