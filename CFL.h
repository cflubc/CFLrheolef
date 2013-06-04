/*
 * CFL.h
 *
 *  Created on: 2013-04-14
 *      Author: ali
 */

#ifndef CFL_H_
#define CFL_H_

#include <string>
#include <sstream>
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


inline void
operator>>( std::istringstream& is, rheolef::point& p )
{p.get(is);}

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

inline void
write_to_diststream( rheolef::odiststream& o )
{}

template< typename T, typename... Args >
inline void
write_to_diststream( rheolef::odiststream& o, std::string const& mark, T const& t, Args const&... args )
{
	o << rheolef::catchmark(mark);
	o << t;
	o << '\n';
	write_to_diststream(o,args...);
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


template< typename Field >
rheolef::Float vector_dot(Field const& f1, Field const& f2 )
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

