/*
 * IncompressibleStokesSolver.h
 *
 *  Created on: 2013-04-13
 *      Author: ali
 */

#ifndef INCOMPRESSIBLESTOKESSOLVER_H_
#define INCOMPRESSIBLESTOKESSOLVER_H_

#include <string>
#include "rheolef.h"
#include "rheolef/diststream.h"

#include "ConfigXML.h"
#include "FlowFields.h"
#include "adaptationCriterions.h"
#include "DiffusionForms.h"

/**
 * Solve incompressible Stokes flow of the form:
 * @f[ -D\Delta u = -\nabla p + f @f]
 * where the @f$ f@f$ term is any general right hand side term provided by user.
 * @f$ D=1@f$ by defualt.
 *
 * The only available BC of this solver is dirichlet for velocity. So may use it
 * for problems like cavity, or problems where the normal stress is zero for all
 * non-dirichelt boudaries of velocity . Say a Poisseulle flow where you apply
 * pressure drop by a constant body force.
 *
 * If you want to have normal stress contribution as well, see
 * NeumannIncompNewtonianStokesSolver class.
 *
 */
template< typename LinearSol >
class IncompLinearDiffusionStokesSolver
{
	typedef rheolef::field field;
	typedef rheolef::geo     geo;

	field& uh;
	field& ph;
	LinearSol solver;

public:
	template< typename FieldsPool, typename DirichletBC >
	IncompLinearDiffusionStokesSolver( const XMLConfigFile& conf,
			                     	   FieldsPool& fields,
			                     	   DirichletBC& BC,
			                     	  const rheolef::Float D=1. ):
		IncompLinearDiffusionStokesSolver( conf, fields,
				     D*rheolef::form(fields.Uh.get_space(),fields.Uh.get_space(),"2D_D") )
	{}

	template< typename FieldsPool, typename DiffusionForm >
	IncompLinearDiffusionStokesSolver( const XMLConfigFile& conf,
			                     	   FieldsPool& fields,
			                     	   const DiffusionForm& Dform
			                     	    ):
		uh(fields.Uh),
		ph(fields.Ph),
		solver( conf, Dform, -rheolef::form(uh.get_space(),ph.get_space(),"div") )
	{}

	void run(){
		solve();
		write_results();
	}

	field adapt_criteria() const
	{return stokes_criterion(uh);}

	void write_results( rheolef::odiststream& o ) const {
		write_field(uh,"u",o);
		write_field(ph,"p",o);
	}

	void write_results() const {
		rheolef::odiststream o(uh.get_geo().name(),"field");
		write_results(o);
		o.close();
	}

	void solve()
	{solver.solve(uh,ph);}

	void solve( const field& rhs )
	{solver.solve(uh,ph,rhs);}

};


#endif /* INCOMPRESSIBLESTOKESSOLVER_H_ */
