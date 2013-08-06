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
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <limits>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "OutputFormatting.h"
#include "TensorFieldIterator.h"
#include "DiffusionForms.h"
#include "ParameterOnFieldDOFs.h"

namespace AugmentedLagrangian_param {

using rheolef::Float;
using rheolef::field;

template< typename Param, typename Initializer >
struct AL_param
{
	typedef Param param_t;
	typedef Initializer init;
};

template< typename Initialzer >
using unique = AL_param<const ParameterOnDOFs::ConstUnique,Initialzer>;

template< typename Initializer >
using non_unique = AL_param<const ParameterOnDOFs::ConstNonUnique,Initializer>;


struct alpha_unique_init
{
	alpha_unique_init( Float const& x ): a(x) {}
	Float get( XMLConfigFile const& conf ) const
	{return a/(a+1.);}
	Float const a;
};

struct alpha_multiRegion_init
{
	alpha_multiRegion_init( Float const& x ): a(x) {}
	void get( field& f, XMLConfigFile const& conf ) const
	{
		typedef ParameterOnDOFs::multiRegion_field_initializer I;
		I::names_t names;
		I::vals_t viscosities;
		I::get_domains_and_values(conf,names,viscosities);
		for(auto& v : viscosities)
			v = a/(v+a);
		I::set_field_regions_values(f,names,viscosities);
	}
	Float const a;
};


typedef unique<alpha_unique_init> alpha_unique;
typedef unique<ParameterOnDOFs::unique_default_initializer> Bingham_unique;

typedef non_unique<alpha_multiRegion_init> alpha_multiRegion;
typedef non_unique<ParameterOnDOFs::multiRegion_field_initializer> Bingham_multiRegion;

} // namespace AugmentedLagrangian_param

template< typename VelocityMinimizationSolver,
          typename BinghamParam,
          typename alphaParam >
class AugmentedLagrangian_basic
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;
	typedef rheolef::space space;
	typedef TensorFieldIterator<field::iterator>       Tensor_itr;
	typedef TensorFieldIterator<field::const_iterator> Tensor_citr;
	typedef Tensor_itr::tensor                         Tensor;

public:
	enum : bool { isVelocityOptimizerLinear = VelocityMinimizationSolver::isLinear };

	template< typename FieldsPool, typename DirichletBC >
	AugmentedLagrangian_basic( const XMLConfigFile& conf,
							   FieldsPool& fields,
							   DirichletBC& BC ):
		Xh(fields.get_geo(), derivative_approx(fields.Uh().get_approx()), "tensor"),
		XML_INIT_VAR(conf,a,"Augmentation_coef"),
		Bn( Xh[0], conf.child({"PhysicalParameters","Bn"}), typename BinghamParam::init() ),
		alpha( Xh[0], conf.child({"PhysicalParameters","Viscosity"}), typename alphaParam::init(a) ),
//		alpha( a/(1.+a) ),
		velocity_minimizer(conf,fields,BC,a),
		Tau(Xh, 0.),
		Gam(Xh, -1000.),
		Gamdot(Xh, 0.),
		TminusaG(Xh, 0.),
		vel_rhs(fields.Uspace(), 0.),
		vel_rhs_const(fields.Uspace(), 0.),
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

	Float iterate_report_strain_change(){
		deltaTau.save_field();
		iterate();
		return deltaTau.calculate_field_change()/a;
	}

	void iterate_ntimes( const int niter){
		for( int i=0; i<niter; ++i )
			iterate();
	}

	void iterate_ntimes_report_strain_velocity_change( const int niter, Float& dGam, Float& dU ){
		iterate_ntimes(niter);
		iterate_report_strain_velocity_change(dGam,dU);
	}

	void iterate_report_strain_velocity_change( Float& dT, Float& dU ){
		save_strain_velocity();
		iterate();
		report_strain_velocity_change(dT,dU);
	}

	void iterate(){
		update_lagrangeMultipliers_fast();
		build_complete_rhs_and_solve_vel_minimization( augmented_lagraniang_rhs() );
	}

	void save_strain_velocity(){
		deltaTau.save_field();
		deltaU.save_field();
	}

	void report_strain_velocity_change( Float& dGam, Float& dU ){
		dGam = deltaTau.calculate_field_change()/a;
		dU = deltaU.calculate_field_change();
	}

	field const augmented_lagraniang_rhs() const
	{return div_ThUh*TminusaG;}

	field adapt_criteria() const {
		// just scalar version of Th
		space T0h( Xh.get_geo(), Xh.get_approx() );
		return interpolate( T0h, sqrt(.5*norm2(Gamdot)+Bn.get_parameter()*sqrt(.5)*norm(Gamdot)) );
	}

	void write_results() const {
		rheolef::odiststream o(Xh.get_geo().name(),"field");
		write_fields(o);
		o.close();
//		o.open("U2StrainRate","m",false);
//		o << rheolef::matlab;
//		Gamdot_server.U2StrainRate.put(o);
//		o.close();
	}

	void write_fields( rheolef::odiststream& o ) const {
		velocity_minimizer.write_results(o);
		write_field(Tau,"T",o);
		write_field(Gam,"Gam",o);
		write_field(Gamdot,"Gammadot",o);
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

	void build_complete_rhs_and_solve_vel_minimization( field const& f ){
		vel_rhs = vel_rhs_const + f;
		solve_vel_minization();
	}

	void solve_vel_minization() const {
		solve_vel_minization(vel_rhs);
	}

	void solve_vel_minization( field const& rhs ) const
	{velocity_minimizer.solve(rhs);}

	field& vel_rhs_const_part()
	{return vel_rhs_const;}

	field& vel_rhs_var_part()
	{return vel_rhs;}

	field const& get_strainRate_lagrangeMultiplier() const
	{return Gam;}

private:
	space Xh;    ///< tensor space of Lagrange multipliers
	Float const a;       ///< Augmentation coef
	typename BinghamParam::param_t Bn;
	typename alphaParam::param_t alpha;
//	Float const alpha;   ///< helper const coeficient

	VelocityMinimizationSolver velocity_minimizer;
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



template< typename VelocityMinimizationSolver, typename BinghamParam, typename alphaParam >
template< typename LoopManipulator >
void AugmentedLagrangian_basic<VelocityMinimizationSolver,BinghamParam,alphaParam>::
update_lagrangeMultipliers( LoopManipulator& obj )
{
	Gamdot_server.get(Gamdot);

	auto B = Bn.begin_dof();
	auto _alpha = alpha.begin_dof();
	Tensor_citr Gdot(Gamdot);
	for( Tensor_itr T(Tau), TmG(TminusaG); !T.end_reached(); ++T, ++TmG, ++Gdot, ++B, ++_alpha, ++obj )
	{
		Tensor TaGdot;
		Float TaGdot_norm(0.);
		for( int i=0; i<Tensor_itr::Ncomp; ++i ){
			TaGdot[i] = T(i)+a*Gdot(i);
			TaGdot_norm += Tensor_itr::coef_for_norm_calc[i]*rheolef::sqr(TaGdot[i]);
		}
		TaGdot_norm = std::sqrt( .5*TaGdot_norm );

		const Float resi = TaGdot_norm-*B;
		if( 0.<resi )
		{
			const Float resi_frac( resi/TaGdot_norm );
			const Float coef = (*_alpha)*resi_frac;
			// update of Lagrange multipliers
			for(int i=0; i<Tensor_itr::Ncomp; ++i){
				T(i)   = TaGdot[i]*(1.-coef);
				TmG(i) = TaGdot[i]*(1.-2.*coef);
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
	assert( TmG.end_reached() );
	assert( Gdot.end_reached() );
}





#endif /* AUGMENTEDLAGRANGIAN_BASIC_H_ */
