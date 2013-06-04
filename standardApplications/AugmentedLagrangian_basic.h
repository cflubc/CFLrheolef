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
#include <string>
#include <iostream>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "OutputFormatting.h"
#include "TensorFieldIterator.h"
#include "DiffusionForms.h"


template< typename VelocityMinimizationSolver >
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
							   DirichletBC& BC ):
		Bn( conf.atof("Bn") ),
		a( conf.atof("a") ),
		alpha( a/(1.+a) ),
		velocity_minimizer(conf,fields,BC,a),
		Xh(fields.get_geo(), derivative_approx(fields.Uh().get_approx()), "tensor"),
		Tau(Xh, 0.),
		Gam(Xh, -1000.),
		Gamdot(Xh, 0.),
		TminusaG(Xh, 0.),
		vel_rhs(fields.Uspace(), 0.),
		vel_rhs_const(vel_rhs.get_space(), 0.),
		Gamdot_server(fields.Uh()),
		div_ThUh( -.5*trans(Gamdot_server.set_desired_strainrate_space(Xh)) ),
		deltaTau(Tau),
		deltaU(fields.Uh())
	{}


	void update_lagrangeMultipliers_fast(){
		const manip_void x;
		update_lagrangeMultipliers(x);
	}

	void update_lagrangeMultipliers_clac_strain_rate_multiplier(){
		StrainRate_lagMultiplier_calc x(Gam,a);
		update_lagrangeMultipliers(x);
	}

	Float iterate_report_stress_change(){
		deltaTau.save_field();
		iterate();
		return deltaTau.calculate_field_change();
	}

	void save_stress_velocity(){
		deltaTau.save_field();
		deltaU.save_field();
	}

	void report_stress_velocity_change( Float& dT, Float& dU ){
		dT = deltaTau.calculate_field_change();
		dU = deltaU.calculate_field_change();
	}

	field const augmented_lagraniang_rhs() const
	{return div_ThUh*TminusaG;}

	field adapt_criteria() const {
		// just scalar version of Th
		space T0h( Xh.get_geo(), Xh.get_approx() );
		return interpolate( T0h, sqrt(.5*norm2(Gamdot)+Bn*std::sqrt(.5)*norm(Gamdot)) );
	}

	void write_results() const {
		rheolef::odiststream o(Xh.get_geo().name(),"field");
		write_fields(o);
		o.close();
	}

	void write_fields( rheolef::odiststream& o ) const {
		velocity_minimizer.write_results(o);
		write_field(Tau,"T",o);
		write_field(Gam,"Gam",o);
	}

	std::string geo_name() const
	{return Xh.get_geo().name();}

	void set_rhs_const_part_to_discrete_dirichlet_rhs(){
		velocity_minimizer.set_discrete_dirichlet_rhs(vel_rhs_const);
	}

	void reset_lagrangeMultipliers(){
		Tau = 0.;
		TminusaG = 0.;
	}

	void iterate_ntimes( const int niter ){
		for( int i=0; i<niter; ++i ){
			update_lagrangeMultipliers_fast();
			build_complete_rhs_and_solve_vel_minimization( augmented_lagraniang_rhs() );
		}
	}

	void iterate_ntimes_report_stress_velocity_change( const int niter, Float& dT, Float& dU ){
		iterate_ntimes(niter);
		iterate_report_stress_velocity_change(dT,dU);
	}

	void iterate_report_stress_velocity_change( Float& dT, Float& dU ){
		save_stress_velocity();
		iterate();
		report_stress_velocity_change(dT,dU);
	}

	void iterate()
	{iterate_ntimes(1);}

	void build_complete_rhs_and_solve_vel_minimization( const field& f ){
		vel_rhs = vel_rhs_const + f;
		solve_vel_minization();
	}

	void solve_vel_minization()
	{velocity_minimizer.solve(vel_rhs);}

	field& vel_rhs_const_part()
	{return vel_rhs_const;}

	field& vel_rhs_var_part()
	{return vel_rhs;}

private:
	Float Bn;      ///< Bingham number
	Float a;       ///< Augmentation parameter
	Float alpha;   ///< helper const coeficient

	VelocityMinimizationSolver velocity_minimizer;
	space Xh;    ///< tensor space of Lagrange multipliers
	field Tau;   ///< Stress Lagrange multiplier
	field Gam;   ///< Strain rate Lagrange multiplier
	field Gamdot;  ///< Strain rate of velocity
	field TminusaG;
	field vel_rhs;
	field vel_rhs_const;
	StrainRateCalculator Gamdot_server;
	rheolef::form div_ThUh;
	L2norm_calculator deltaTau;
	L2norm_calculator deltaU;

	template< typename LoopManipulator >
	void update_lagrangeMultipliers( LoopManipulator& obj );

	struct manip_void {
		void operator++() const {}
		void do_extra_stuff( int, Float, Float ) const {}
	};

	struct StrainRate_lagMultiplier_calc {
		StrainRate_lagMultiplier_calc( field& Gamma, Float const& _a ):
			G(Gamma),
			beta(1./(1.+_a))
		{}
		void operator++() {++G;}
		void do_extra_stuff( int const i, Float const& TaG, Float const& resi_stress_frac )
		{G(i) = TaG*resi_stress_frac*beta;}
		Tensor_itr G;
		Float const beta;
	};
};



template< typename VelocityMinimizationSolver>
template< typename LoopManipulator >
void AugmentedLagrangian_basic<VelocityMinimizationSolver>::
update_lagrangeMultipliers( LoopManipulator& obj )
{
	assert_equal(Gamdot,Tau);
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
			TaGdot_norm += Tensor_itr::coef_for_norm_calc[i]*rheolef::sqr(TaGdot[i]);
		}
		TaGdot_norm = std::sqrt( .5*TaGdot_norm );

		const Float resi = TaGdot_norm-Bn_;
		if( 0.<resi )
		{
			const Float resi_frac( resi/TaGdot_norm );
			const Float coef = alpha_*resi_frac;
			// update of Lagrange multipliers
			for(int i=0; i<Tensor_itr::Ncomp; ++i){
				T(i)   = TaGdot[i]*(1.-coef);
				TmG(i) = TaGdot[i]*(1.-coef-coef);
				obj.do_extra_stuff(i,TaGdot[i],resi_frac);
			}
		}
		else
		{
			// unyielded, update of Lagrange multipliers
			for(int i=0; i<Tensor_itr::Ncomp; ++i){
				T(i)   = TaGdot[i];
				TmG(i) = TaGdot[i];
				obj.do_extra_stuff(i,TaGdot[i],0.);
			}
		}
	}
}





#endif /* AUGMENTEDLAGRANGIAN_BASIC_H_ */
