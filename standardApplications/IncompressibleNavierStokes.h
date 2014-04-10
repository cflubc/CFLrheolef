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

#include "rheolef/diststream.h"

#include "IncompressibleNavierStokes_core.h"
#include "CFL.h"
#include "ConfigXML.h"
#include "ResidualTablePrinter.h"
#include "GenericIteration.h"
#include "ErrorAnalysis.h"


class IncompressibleNavierStokes : public IncompressibleNavierStokes_core
{
public:

	template< typename FieldsPool, typename DirichletBC >
	IncompressibleNavierStokes( XMLConfigFile const& conf, FieldsPool& fields, DirichletBC& BC, Float const viscosity=1. ):
		IncompressibleNavierStokes_core(conf,fields,BC,viscosity),
		deltaU(uh1,uh2),
		output(30,std::cout,10,16),
		converge_loop(conf)
	{}

	void run()
	{
		solver.set_discrete_dirichlet_rhs(dirichlet_rhs,uh);
		uh1 = uh;
		converge_loop(*this);
		rheolef::odiststream o (uh.get_geo().name(), "field");
		write_to_diststream(o,"u",uh, "p",ph);
		o.close();
	}

	void iterate()
	{
		compute_convective_rhs();
		field const lh = dirichlet_rhs + convect_rhs;
		solve(lh);
	}

	void iterate_report( size_t const niter, Float& res )
	{
		deltaU.save_field();
		iterate();
		res = deltaU.calculate_field_change()/delta_t;
		output.print_header_if_needed("iteration","|Un+1-Un/dt|L2");
		output.print(niter,res);
	}

	field adapt_criteria() const
	{
		space T0h( uh.get_geo(), derivative_approx(uh.get_approx()) );
		return interpolate( T0h, sqrt(Re*norm2(uh) + 4.*norm2(D(uh))) );
	}

	L2norm_calculator deltaU;
	ResidualTablePrinter<2,std::ostream> output;
	GenericIteration converge_loop;
};

#endif /* INCOMPRESSIBLENAVIERSTOKES_H_ */
