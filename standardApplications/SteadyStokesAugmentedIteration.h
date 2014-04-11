/*
 * SteadyStokesAugmentedIteration.h
 *
 *  Created on: Sep 15, 2013
 *      Author: ali
 */

#ifndef STEADYSTOKESAUGMENTEDITERATION_H_
#define STEADYSTOKESAUGMENTEDITERATION_H_

#include <cstddef>
#include <iostream>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "ErrorAnalysis.h"
#include "GenericIteration.h"
#include "ResidualTablePrinter.h"
#include "BlockSystem_abtb.h"
#include "IncompressibleStokesSolver.h"


template< typename BasicAugmentedLagrangian >
class SteadyStokesAugmentedIteration
{
	typedef std::size_t size_t;
	typedef rheolef::Float Float;
	typedef rheolef::field field;

public:
	typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> VelocityMinimizer;

	template< typename FieldsPool, typename DirichletBC >
	SteadyStokesAugmentedIteration( const XMLConfigFile& conf, FieldsPool& fields, DirichletBC& BC ):
		AL(conf,fields,BC),
		velocity_minimizer(conf,fields,BC,AL.augmentation_coef()),
		deltaU(fields.Uh()),
		algo(conf),
		vel_rhs_const_part(fields.Uspace(), 0.),
		vel_rhs(vel_rhs_const_part),
		base_name( fields.geo_name() ),
		output(30,std::cout,10,16,16)
	{}

	void write_results()
	{
		rheolef::odiststream o(base_name,"field");
		velocity_minimizer.write_results(o);
		AL.write_fields(o);
		o.close();
	}

	void do_iterations( field const& physical_rhs ){
		velocity_minimizer.set_discrete_dirichlet_rhs(vel_rhs_const_part);
		vel_rhs_const_part += physical_rhs;
		algo(*this);
	}

	field adapt_criteria() const
	{return AL.adapt_criteria();}

	void iterate(){
		solve_vel();
		AL.update_lagrangeMultipliers_fast();
	}

	void iterate_report( size_t const niter, Float& res)
	{
		deltaU.save_field();
		solve_vel();
		Float const Gamres = AL.update_lagrangeMultipliers_report_strain_residual();
		Float const Ures = deltaU.calculate_field_change();
		res = rheolef::max(Gamres,Ures);
		output.print_header_if_needed("\niteration","|Un-Un-1|L2","|Gam-Gamdot|L2");
		output.print(niter,Ures,Gamres);
	}

	field get_strainRate_lagrangeMultiplier() const
	{return AL.get_strainRate_lagrangeMultiplier();}

private:

	void solve_vel(){
		vel_rhs = vel_rhs_const_part + AL.augmented_lagraniang_rhs();
		velocity_minimizer.solve(vel_rhs);
	}

	BasicAugmentedLagrangian AL;
	VelocityMinimizer velocity_minimizer;
	L2norm_calculator deltaU;
	GenericIteration algo;
	field vel_rhs_const_part;
	field vel_rhs;
	std::string const base_name;
	ResidualTablePrinter<3,std::ostream> output;
};



template< typename BasicAugmentedLagrangian, typename VelocityRHSManipulator >
class SteadyStokesAugmentedLagrangian
{
	SteadyStokesAugmentedIteration<BasicAugmentedLagrangian> iter;
	VelocityRHSManipulator physical_rhs_term;

public:

	template< typename FieldsPool, typename DirichletBC >
	SteadyStokesAugmentedLagrangian( const XMLConfigFile& conf, FieldsPool& fields, DirichletBC& BC ):
		iter(conf, fields, BC),
		physical_rhs_term( conf.child("source_term"), fields.Uspace() )
	{}

	void run()
	{
		iter.do_iterations( physical_rhs_term.get_rhs() );
		iter.write_results();
	}

};

#endif /* STEADYSTOKESAUGMENTEDITERATION_H_ */
