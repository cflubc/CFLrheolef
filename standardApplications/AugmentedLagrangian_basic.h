/*
 * AugmentedLagrangian.h
 *
 *  Created on: 2013-04-19
 *      Author: ali
 */


//	//no idea why when use const& for return type I get segfault...
//	static const XMLConfigFile ALconf( const XMLConfigFile& conf )
//	{return conf.child("AugmentedLagrangian"); }

#ifndef AUGMENTEDLAGRANGIAN_BASIC_H_
#define AUGMENTEDLAGRANGIAN_BASIC_H_


#include <cmath>
#include <cstdio>
#include <iostream>
#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "OutputFormatting.h"
#include "TensorFieldIterator.h"
#include "DiffusionForms.h"


template< typename FlowSolver >
class AugmentedLagrangian_basic
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;
	typedef rheolef::space space;
	typedef TensorFieldIterator<field::iterator>       Tensor_itr;
	typedef TensorFieldIterator<field::const_iterator> Tensor_citr;
	typedef Tensor_itr::tensor                         Tensor;

public:
	template< typename FieldsPool, typename DirichletBC >
	AugmentedLagrangian_basic( const XMLConfigFile& conf,
							   FieldsPool& fields,
							   DirichletBC& BC
							   ):
		Bn( conf.atof("Bn") ),
		a( conf.atof("a") ),
		alpha( a/(1.+a) ),
		flow_solver(conf,fields,BC,a),
		Xh(fields.Uh.get_geo(), derivative_approx(fields.Uh.get_approx()), "tensor"),
		Tau(Xh, 0.),
		Gam(Xh, 0.),
		Gamdot(Xh, 0.),
		TminusaG(Xh, 0.),
		Gamdot_server(fields.Uh),
		div_ThUh( -.5*trans(Gamdot_server.set_desired_strainrate_space(Xh)) ),
		deltaTau(Tau)
	{}


	void contribute_to_rhs_fast( field& rhs ){
		const manip_void x;
		add_rhs(rhs,x);
	}

	Float contribute_to_rhs_report_stress_change( field& rhs ){
		deltaTau.save_field();
		contribute_to_rhs_fast(rhs);
		return deltaTau.calculate_field_change();
	}

	field adapt_criteria() const {
		// just scalar version of Th
		space T0h( Xh.get_geo(), Xh.get_approx() );
		return interpolate( T0h, sqrt(norm2(Gamdot)+Bn*norm(Gamdot)) );
	}

	void write_results() const {
		rheolef::odiststream o(Xh.get_geo().name(),"field");
		write_fields(o);
		o.close();
	}

	void write_fields( rheolef::odiststream& o ) const {
		flow_solver.write_results(o);
		write_field(Tau,"T",o);
		write_field(Gam,"Gam",o);
	}

	void solve( field& rhs ){
		flow_solver.solve(rhs);
		rhs = 0.;
	}


private:
	Float Bn;      ///< Bingham number
	Float a;       ///< Augmentation parameter
	Float alpha;   ///< helper variable

	FlowSolver flow_solver;
	space Xh;    ///< tensor space of Lagrange multipliers
	field Tau;   ///< Stress Lagrange multiplier
	field Gam;   ///< Strain rate Lagrange multiplier
	field Gamdot;  ///< Strain rate of velocity
	field TminusaG;
	StrainRateCalculator Gamdot_server;
	rheolef::form div_ThUh;
	L2norm_calculator deltaTau;

	template< typename LoopManipulator >
	void add_rhs( field& rhs, LoopManipulator& obj );

	struct manip_void {
		void operator++() const {}
		void set_stress_residual_fraction( Float ) const {}
		void do_extra_stuff( int, Float, Float ) const {}
	};
};



template< typename FlowSolver>
template< typename LoopManipulator >
void AugmentedLagrangian_basic<FlowSolver>::add_rhs( field& rhs, LoopManipulator& obj )
{
	assert_equal(Gamdot,Tau);
	//Tensor_itr G(Gam); //for last iteration when want to save...
	const Float& alpha_ = alpha;
	const Float& Bn_ = Bn;
	const Float& a_ = a;

	Gamdot_server.get(Gamdot);
	Tensor_citr Gdot(Gamdot);
	for( Tensor_itr T(Tau), TmG(TminusaG); !T.end_reached(); ++T, ++TmG, ++Gdot, ++obj )
	{
		Tensor TaGdot;
		Float TaGdot_norm(0.);
		for( int i=0; i<Tensor_itr::Ncomp; ++i ){
			TaGdot[i] = T(i)+a_*Gdot(i);
			TaGdot_norm += Tensor_itr::coef_for_norm_calc[i]*rheolef::sqr(TaGdot[i]); //*TaGdot[i];
		}
		TaGdot_norm = std::sqrt( .5*TaGdot_norm );

		const Float resi = TaGdot_norm-Bn_;
		Float coef(0.);
		if( 0.<resi ){
			const Float& resi_frac( resi/TaGdot_norm );
			obj.set_stress_residual_fraction(resi_frac);
			coef = alpha_*resi_frac;
		}

		// update of Lagrange mults
		for(int i=0; i<Tensor_itr::Ncomp; ++i){
			T(i)   = TaGdot[i]*(1.-coef);
			TmG(i) = TaGdot[i]*(1.-2.*coef);
			obj.do_extra_stuff( i, TaGdot[i], T(i) );
//				G(i)   = TaGdot[i]*coef;
		}
	}
	rhs += div_ThUh*TminusaG;
}





#endif /* AUGMENTEDLAGRANGIAN_BASIC_H_ */
