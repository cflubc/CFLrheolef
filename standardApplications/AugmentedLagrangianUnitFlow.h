/*
 * AugmentedLagrangianUnitFlow.h
 *
 *  Created on: 2013-05-17
 *      Author: ali
 */

#ifndef AUGMENTEDLAGRANGIANUNITFLOW_H_
#define AUGMENTEDLAGRANGIANUNITFLOW_H_

#include <cstdlib>
#include <cmath>
#include "rheolef.h"

#include "ConfigXML.h"
#include "ErrorAnalysis.h"
#include "AugmentedLagrangian_basic.h"

template< typename T = double >
class SecantMethod
{
	const T target;
	T f1;
	T x2;
	T dx;
	T tolerance;
	size_t max_iter;
	size_t iter;
	bool converged;

public:
	T predict_new_input( T const& f2 )
	{
		T const dtarget = target-f2;
		converged = rheolef::abs(dtarget)<tolerance;
		if( !converged ){
			++iter;
			dx = dx/(f2-f1)*dtarget;
			// update to new x
			x2 += dx;
			f1 = f2;
		}
		return x2;
	}

	SecantMethod( size_t nmax, T const& tol, T const& _target, T const& x0, T const& f0, T const& x ):
		target(_target),
		f1(f0),
		x2(x),
		dx(x-x0),
		tolerance(tol),
		max_iter(nmax),
		iter(0)
	{}

	void reset(){
		iter = 0;
		converged = false;
	}

	size_t n_iterations_done() const
	{return iter;}

	bool not_converged_and_have_iterations_left() const {
		return (iter<max_iter) && !converged;
	}

};


template< typename VelocityMinimizationSolver, typename FlowControler >
class AugmentedLagrangianUnitFlow
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

	field& uh;
public:
	template< typename FieldsPool, typename DirichletBC >
	AugmentedLagrangianUnitFlow( const XMLConfigFile& conf,
								 FieldsPool& fields,
								 DirichletBC BC ):
		flow_control(conf,fields),
		AL(conf,fields,BC),
		Uchange(fields.Uh),
		uh(fields.Uh)
	{}

	template< typename IndirectField >
	static Float vector_dot(IndirectField const& f1, IndirectField const& f2 ){
		Float s = 0.;
		for( auto i1=f1.begin_dof(), i2=f2.begin_dof(); i1!=f1.end_dof(); ++i1,++i2 ){
			s += (*i1)*(*i2);
		}
		return s;
	}

	class BoundaryFluxCalculator
	{
		std::string const boundary_name;
		field const& f;
		field m;

	public:
		BoundaryFluxCalculator( const std::string& name, const field& _f ):
			boundary_name(name),
			f(_f)
		{
			rheolef::test v( f.get_space() );
			m = integrate( boundary_name, dot(rheolef::normal(),v) );
		}

		Float calc_flux() const {
			return vector_dot( m[boundary_name], f[boundary_name] );
		}
	};

	void run()
	{
//		field m = integrate("left",dot(rheolef::normal(),v));
		rheolef::test v{uh.get_space()};
		BoundaryFluxCalculator inlet("left",uh);
		field bodyf = integrate( dot(rheolef::point(1.,0.),v) );
		//flow_controler ctor(conf,uh)
		// contains fluxcalc, predictor

		int niter = 0;
		AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
		//flowrate controller
//		flow_control.get_control_parameter();

		Float cparam {4.};
		SecantMethod<Float> predictor(6,1E-5,2.,0.,0.,cparam);
		do {
			predictor.reset();
			while( predictor.not_converged_and_have_iterations_left() ){
				printf(" %f,",cparam);
				AL.vel_rhs_var_part() = cparam*bodyf; //flow_control.set_control_rhs()
				AL.solve_vel_minization();
				const Float flux = rheolef::abs( inlet.calc_flux() ); //flow_control.calculate_flux()
				cparam = predictor.predict_new_input(flux);
			}

			AL.set_rhs_const_part_to_discrete_dirichlet_rhs();
			const Float Tres = AL.set_rhs_report_stress_change();
			printf("  Tres:%e\n",Tres);
			AL.vel_rhs_const_part() += AL.vel_rhs_var_part();
			++niter;
		} while( niter<max_iteration_stage1 );
		AL.write_results();
	}

private:
	FlowControler flow_control;
	AugmentedLagrangian_basic<VelocityMinimizationSolver> AL;

	int max_iteration_stage1 = 100;
	int iterations_before_unitflow {1};
	L2norm_calculator Uchange;
};


#endif /* AUGMENTEDLAGRANGIANUNITFLOW_H_ */
