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
#include "ConfigXML.h"


class IncompressibleNavierStokes
{
	typedef rheolef::Float Float;
	typedef rheolef::field field;
	typedef rheolef::space space;
	typedef rheolef::form form;
	typedef rheolef::characteristic characteristic;
	typedef rheolef::odiststream odiststream;
	typedef std::string string;

	int navier_stokes_solve (field l0h, odiststream *p_derr=0)
	{
	  const space& Xh = uh.get_space();
	  const space& Qh = ph.get_space();
	  string label = "navier-stokes-" + Xh.get_geo().name();
	  rheolef::quadrature_option_type qopt;
	  qopt.set_family(rheolef::quadrature_option_type::gauss_lobatto);
	  qopt.set_order(Xh.degree());
	  rheolef::trial u (Xh), p (Qh);
	  rheolef::test  v (Xh), q (Qh);
	  form mp = integrate (p*q, qopt);
	  form m  = integrate (dot(u,v), qopt);
	  form a  = integrate (2*ddot(D(u),D(v)) + 1.5*(Re/delta_t)*dot(u,v), qopt);
	  form b  = integrate (-div(u)*q, qopt);
	  rheolef::solver_abtb stokes (a.uu(), b.uu(), mp.uu());
	  if (p_derr != 0) *p_derr << "[" << label << "] #n |du/dt|" << std::endl;
	  field uh1 = uh;
	  for (size_t n = 0; true; n++) {
	    field uh2 = uh1;
	    uh1  = uh;
	    field uh_star = 2.0*uh1 - uh2;
	    characteristic X1 (    -delta_t*uh_star);
	    characteristic X2 (-2.0*delta_t*uh_star);
	    field l1h = integrate (dot(compose(uh1,X1),v), qopt);
	    field l2h = integrate (dot(compose(uh2,X2),v), qopt);
	    field lh  = l0h + (Re/delta_t)*(2*l1h - 0.5*l2h);
	    stokes.solve (lh.u() - a.ub()*uh.b(), -(b.ub()*uh.b()),
	                  uh.set_u(), ph.set_u());


	    field duh_dt = (3*uh - 4*uh1 + uh2)/(2*delta_t);
	    Float residual = sqrt(m(duh_dt,duh_dt));
	    if (p_derr != 0) *p_derr << "[" << label << "] "<< n << " " << residual << std::endl;
	    if (residual < tol) {
	      tol = residual;
	      max_iter = n;
	      return 0;
	    }
	    if (n == max_iter-1) {
	      tol = residual;
	      return 1;
	    }
	  }
	}
public:
	enum : bool { isLinear=false };

	template< typename FieldsPool, typename DirichletBC >
	IncompressibleNavierStokes( XMLConfigFile const& conf, FieldsPool& fields, DirichletBC& BC ):
	uh(fields.Uh()),
	ph(fields.Ph()),
	Re( conf({"PhysicalParameters","Re"},Re) ),
	delta_t( conf({"PhysicalParameters","dt"},delta_t) ),
	XML_INIT_VAR(conf,tol,"tolerance"),
	XML_INIT_VAR(conf,max_iter,"max_iteration")
	{}

	void run()
	{
		using rheolef::catchmark;
		field fh(uh.get_space(),0);
		navier_stokes_solve (fh, &rheolef::derr);
		odiststream o (uh.get_geo().name(), "field");
		o << catchmark("Re") << Re << std::endl
		  << catchmark("delta_t") << delta_t << std::endl
		  << catchmark("u")  << uh
		  << catchmark("p")  << ph;
		o.close();
	}

	field adapt_criteria() const
	{}

	Float const Re;
	Float const delta_t;
	Float tol;
	std::size_t max_iter;

	field& uh;
	field& ph;
};


#endif /* INCOMPRESSIBLENAVIERSTOKES_H_ */
