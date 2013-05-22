/*
 * CFL.h
 *
 *  Created on: 2013-04-14
 *      Author: ali
 */

#ifndef CFL_H_
#define CFL_H_

#include <string>
#include <initializer_list>

#include "rheolef.h"
#include "rheolef/diststream.h"


typedef const char* cstr;

/**
 * The type used to show that the function is getting several objects
 */
template< typename T >
using variadic_input = std::initializer_list< const T >;


inline void
write_field( const rheolef::field& f, const char* mark, rheolef::odiststream& o )
{
	o << rheolef::catchmark(mark);
	o << f;
}

void assert_equal( const rheolef::field& f1, const rheolef::field& f2 );
std::string derivative_approx( const std::string& approx );
void print_solution_convergence_message( bool converged );


struct StrainRateCalculator
{
	const rheolef::field& Uh;
	rheolef::form U2StrainRate;

	StrainRateCalculator():
		Uh(rheolef::field())
	{}

	StrainRateCalculator( const rheolef::field& u ):
		Uh(u)
	{}

	rheolef::form set_desired_strainrate_space( const rheolef::space& Th )
	{
		rheolef::form _2D( Uh.get_space(), Th, "2D" );
		rheolef::form inv_mt( Th, Th, "inv_mass" );
		U2StrainRate = inv_mt*_2D;
		return _2D;
	}

	void get( rheolef::field& strain_rate )
	{ strain_rate = U2StrainRate*Uh; }
};


class RecuringAlarm
{
	int n_recuring;
	int counter;

public:
	RecuringAlarm( int n, int initval=0 ):
		n_recuring(n),
		counter(initval)
	{}

	bool alarm_ringing()
	{ return (counter++%n_recuring)==0; }

	void reset()
	{ counter=0; }
};


template< typename IndirectField >
rheolef::Float vector_dot(IndirectField const& f1, IndirectField const& f2 )
{
	rheolef::Float s = 0.;
	for( auto i1=f1.begin_dof(), i2=f2.begin_dof(); i1!=f1.end_dof(); ++i1,++i2 )
		s += (*i1)*(*i2);
	return s;
}


#endif /* CFL_H_ */

