/*
 * SteadyNavierStokesAugmentedLagrangian.h
 *
 *  Created on: Oct 8, 2013
 *      Author: ali
 */

#ifndef STEADYNAVIERSTOKESAUGMENTEDLAGRANGIAN_H_
#define STEADYNAVIERSTOKESAUGMENTEDLAGRANGIAN_H_

#include <iostream>

#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "GenericIteration.h"
#include "ResidualTablePrinter.h"
#include "IncompressibleNavierStokes_core.h"


template< typename BasicAugmentedLagrangian, typename VelocityRHSManipulator  >
class SteadyNavierStokesAugmentedLagrangian : public IncompressibleNavierStokes_core
{
public:

	template< typename FieldsPool, typename DirichletBC >
	SteadyNavierStokesAugmentedLagrangian( const XMLConfigFile& conf, FieldsPool& fields, DirichletBC& BC ):
		IncompressibleNavierStokes_core(conf,fields,BC,conf("Augmentation_coef",Float())),
		AL(conf,fields,BC),
		physical_rhs_term( conf.child("source_term"), fields.Uspace() ),
		deltaU(uh1,uh2),
		converge_loop(conf),
		vel_rhs_const_part(fields.Uspace(), 0.),
		vel_rhs(vel_rhs_const_part),
		output(30,std::cout,10,20,20)
	{}

	void run()
	{
		solver.set_discrete_dirichlet_rhs(dirichlet_rhs,uh);
		vel_rhs_const_part = dirichlet_rhs; // + physical_rhs_term.get_rhs();
		uh1 = uh;
		converge_loop(*this);
		AL.update_lagrangeMultipliers_clac_strain_rate_multiplier();
		rheolef::odiststream o (uh.get_geo().name(), "field");
		write(o);
		AL.write_fields(o);
		o.close();
	}

	void iterate(){
		compute_convective_rhs();
		AL.update_lagrangeMultipliers_fast();
		vel_rhs = vel_rhs_const_part + convect_rhs + AL.augmented_lagraniang_rhs();
		solve(vel_rhs);
	}

	void iterate_report( size_t const niter, Float& res )
	{
		deltaU.save_field();
		AL.save_strain();
		iterate();
		Float const Gamres = AL.strain_change()/delta_t;
		Float const Ures = deltaU.calculate_field_change()/delta_t;
		res = rheolef::max(Gamres,Ures);
		output.print_header_if_needed("iteration","|Un+1-Un/dt|L2","|Gam-Gamdot/dt|L2");
		output.print(niter,Ures,Gamres);
	}

	field adapt_criteria() const {
		space T0h( uh.get_geo(), derivative_approx(uh.get_approx()) );
		return interpolate( T0h, sqrt(Re*norm2(uh) + norm2(AL.Gamdot)
				                + AL.Bn.get_parameter()*norm(AL.Gamdot)) ); //4.*norm2(D(uh))
	}

	BasicAugmentedLagrangian AL;
	VelocityRHSManipulator physical_rhs_term;
	L2norm_calculator deltaU;
	GenericIteration converge_loop;
	field vel_rhs_const_part;
	field vel_rhs;
	ResidualTablePrinter<3,std::ostream> output;
};


#endif /* STEADYNAVIERSTOKESAUGMENTEDLAGRANGIAN_H_ */
