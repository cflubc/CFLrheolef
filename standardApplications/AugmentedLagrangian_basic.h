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

template< typename BinghamParam, typename alphaParam >
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
		Xh(fields.get_geo(), derivative_approx(fields.Uh().get_approx()), "tensor"),
		XML_INIT_VAR(conf,a,"Augmentation_coef"),
		Bn( Xh[0], conf.child({"PhysicalParameters","Bn"}), typename BinghamParam::init() ),
		alpha( Xh[0], conf.child({"PhysicalParameters","Viscosity"}), typename alphaParam::init(a) ),
		Tau(Xh, 0.),
		Gam(Xh, -1000.),
		Gamdot(Xh, 0.),
		TminusaG(Xh, 0.),
		Gamdot_server(fields.Uh()),
		div_ThUh( -.5*trans(Gamdot_server.set_desired_strainrate_space(Xh)) ),
		deltaTau(Tau)
	{}


	void update_lagrangeMultipliers_fast(){
		const manip_void x;
		update_lagrangeMultipliers(x);
	}

	void update_lagrangeMultipliers_clac_strain_rate_multiplier(){
		StrainRate_lagMultiplier_calc<typename alphaParam::param_t::const_iterator> x(Gam,alpha.begin_dof(),a);
		update_lagrangeMultipliers(x);
	}

	Float update_lagrangeMultipliers_report_strain_residual(){
		save_strain();
		update_lagrangeMultipliers_fast();
		return strain_change();
	}

	void save_strain()
	{deltaTau.save_field();}

	Float strain_change()
	{return deltaTau.calculate_field_change()/a;}

	field augmented_lagraniang_rhs() const
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
		write_field(Tau,"T",o);
		write_field(Gam,"Gam",o);
		write_field(Gamdot,"Gammadot",o);
	}

	std::string geo_name() const
	{return Xh.get_geo().name();}

	void reset_lagrangeMultipliers(){
		Tau = 0.;
		TminusaG = 0.;
	}

	field const& get_strainRate_lagrangeMultiplier() const
	{return Gam;}

	Float augmentation_coef() const
	{return a;}


	space Xh;    ///< tensor space of Lagrange multipliers
	Float const a;       ///< Augmentation coef
	typename BinghamParam::param_t Bn;
	typename alphaParam::param_t alpha;

	field Tau;   ///< Stress Lagrange multiplier
	field Gam;   ///< Strain rate Lagrange multiplier
	field Gamdot;  ///< Strain rate of velocity
	field TminusaG;
	StrainRateCalculator Gamdot_server;
	rheolef::form div_ThUh;
	L2norm_calculator deltaTau;

	template< typename LoopManipulator >
	void update_lagrangeMultipliers( LoopManipulator& obj );

	struct manip_void {
		void operator++() const {}
		void do_extra_stuff( int, Float, Float ) const {}
	};

	template< typename alphaIterator >
	struct StrainRate_lagMultiplier_calc {
		StrainRate_lagMultiplier_calc( field& Gamma, alphaIterator const& itr, Float const& _a ):
			G(Gamma),
			a_coef(_a),
			alpha(itr)
		{}
		void operator++() {
			++G;
			++alpha;
		}
		void do_extra_stuff( int const i, Float const& TaG, Float const& resi_stress_frac )
		{G(i) = TaG*resi_stress_frac*(*alpha)/a_coef;}

		Tensor_itr G;
		Float const& a_coef;
		alphaIterator alpha;
	};
};



template< typename BinghamParam, typename alphaParam >
template< typename LoopManipulator >
void AugmentedLagrangian_basic<BinghamParam,alphaParam>::
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
