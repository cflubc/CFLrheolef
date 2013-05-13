/*
 * ConvergenceMonitor.h
 *
 *  Created on: 2013-04-26
 *      Author: ali
 */

#ifndef CONVERGENCEMONITOR_H_
#define CONVERGENCEMONITOR_H_

#include <cstdlib>
#include <vector>
#include <string>
#include <initializer_list>


class ConvergenceMonitor
{
	typedef std::string string;
	typedef const char* cstr;
	typedef std::vector<double> convergence_history;

	string file_name_base;
	double absolute_error;
	std::vector<size_t> iteration;
	std::vector<cstr> param_names;
	std::vector<convergence_history> converge_histories;

	size_t n_parameters() const
	{ return param_names.size(); }

public:
	ConvergenceMonitor( const string& name,
			            const double& error_limit,
			            std::initializer_list<cstr> names );

	void add_point( const size_t& iter, std::initializer_list<double> vals );
	bool is_converged() const;

	/**
	 * Save convergence history in file. Name of saved file is the given
	 * name + an optional suffix, which can be usefull in saving different
	 * mesh adaptation convergence for the same parameter, say:
	 * U.cvg, U-1.cvg, ...
	 *
	 * This function also generates eps plot of parameters(gnuplot)
	 */
	void save_to_file( const string& suffix="" ) const;
};



#endif /* CONVERGENCEMONITOR_H_ */

