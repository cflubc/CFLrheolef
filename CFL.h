/*
 * CFL.h
 *
 *  Created on: 2013-04-14
 *      Author: ali
 */

#ifndef CFL_H_
#define CFL_H_

#include <string>
#include <iomanip>
#include <initializer_list>

#include "rheolef.h"
#include "rheolef/diststream.h"

#include "MemoryUseage.h"
#include "PrintArguments.h"

typedef const char* cstr;

constexpr cstr CFL_FieldsPool_Module = "FEfields";
constexpr cstr CFL_SaveFolder_BaseName = "result";
constexpr static double PI = std::acos(-1);

inline std::string
domain_filename( std::string const& base )
{return base+".dmn";}

inline std::string
geo_filename( std::string const& base )
{return base+".geo";}

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

void CFL_mkresult_folder_and_cd_to_it( int iadapt );


template< typename Stream >
void CFL_print_time_memory_useage( Stream& s, double const& runtime  )
{
	s << std::setprecision(3);
	s << "\n--------------------------------\n";
	print_args(s,"Running time: ",runtime," hours\n");
	print_memory_useage( s );
	s <<   "--------------------------------\n";
}


class XMLConfigFile;
/**
 * put no or remove <plot_mesh_args> part for not plotting the mesh.
 * all args are sent to rheolef "geo" command. So look at geo manual for
 * possible arguments. Most important ones are:
 *  1. "-gnuplot":  to view mesh interactively
 *  2. "-image-format eps" : to save mesh eps file
 * @param conf
 */
void plot_mesh( XMLConfigFile const& conf, std::string const& geofile );

#endif /* CFL_H_ */

