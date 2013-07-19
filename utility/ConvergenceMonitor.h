/*
 * ConvergenceMonitor.h
 *
 *  Created on: 2013-04-26
 *      Author: ali
 */

#ifndef CONVERGENCEMONITOR_H_
#define CONVERGENCEMONITOR_H_

#include <cstddef>
#include <vector>
#include <string>
#include <initializer_list>


class ConvergenceMonitor
{
	typedef std::string string;
	typedef std::size_t size_t;
	typedef const char* cstr;

	size_t n_parameters() const
	{ return param_names.size(); }

public:

	typedef std::vector<double> convergence_history;

	ConvergenceMonitor( string const& file_name,
			            std::initializer_list<cstr> param_names );

	void add_point( const size_t iter, std::initializer_list<double> vals );

	bool is_converged( double const& err ) const;

	convergence_history const& operator[]( size_t i )
	{return converge_histories[i];}

	void rename_and_init( string const& name, std::initializer_list<cstr> names );
	void rename_and_init( string const& name );

	void clear();

	/**
	 * Save convergence history in file. Name of saved file is the given
	 * name + an optional suffix, which can be usefull in saving different
	 * mesh adaptation convergence for the same parameter, say:
	 * U.cvg, U-1.cvg, ...
	 *
	 * This function also generates eps plot of parameters(gnuplot)
	 */
	void save_to_file( const string& suffix="" ) const;


private:

	string file_name_base;
	std::vector<size_t> iteration;
	std::vector<cstr> param_names;
	std::vector<convergence_history> converge_histories;
};


#endif /* CONVERGENCEMONITOR_H_ */

