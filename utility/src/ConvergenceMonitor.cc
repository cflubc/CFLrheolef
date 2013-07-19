/*
 * ConvergenceMonitor.cc
 *
 *  Created on: 2013-04-27
 *      Author: ali
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>

#include "ConvergenceMonitor.h"


ConvergenceMonitor::ConvergenceMonitor(
		 	 	 	string const& name,
					std::initializer_list<cstr> names ):
	file_name_base(name),
	param_names( begin(names), end(names) ),
	converge_histories( names.size() )
{}


void
ConvergenceMonitor::add_point( const size_t iter, std::initializer_list<double> vals )
{
	if( vals.size()!=converge_histories.size() )
		throw std::logic_error("The number of provided points is different from class initialization");

	iteration.push_back(iter);
	auto val( begin(vals) );
	for(auto& h : converge_histories){
		h.push_back(*val);
		++val;
	}
}


bool
ConvergenceMonitor::is_converged( double const& err ) const
{
	bool ans(true);
	for(const auto& h : converge_histories)
		ans = ans && (h.back()<=err);
	return ans;
}


void
ConvergenceMonitor::rename_and_init( string const& name,
					     			 std::initializer_list<cstr> names )
{
	file_name_base = name;
	param_names.assign( begin(names), end(names) );
	clear();
}

void
ConvergenceMonitor::rename_and_init( string const& name )
{
	file_name_base = name;
	clear();
}

void
ConvergenceMonitor::clear()
{
	iteration.clear();
	converge_histories.clear();
	converge_histories.resize( n_parameters() );
}


void
ConvergenceMonitor::save_to_file( const string& suffix ) const
{
	const string full_name(file_name_base+suffix);

	const string out_file(full_name+".cvg");
	std::ofstream o(out_file);
	// header of file
	o << "# iteration\t";
	for(auto& name : param_names){
		o << name;
		o << '\t';
	}
	o << '\n';
	// writing parameters row by row
	for(size_t i(0); i<iteration.size(); ++i){
		o << iteration[i];
		o << '\t';
		for(const auto& h : converge_histories){
			o << h[i];
			o << '\t';
		}
		o << '\n';
	}
	o.close();

	// generate convergence plot using gnuplot
	const string gnuplot_file(full_name+".gnu");
	o.open(gnuplot_file);
	o << "set term postscript eps\n"
		 "set logscale y\n"
		 "set xlabel 'iteration'\n"
		 "set output '"+full_name+".eps'\n\n";

	o << "plot \\\n";
	for(size_t i=0; i<n_parameters(); ++i){
		cstr tail;
		if( i<n_parameters()-1 )
			tail = ", \\";
		else // last line no need for ,
			tail = " ";
		char line[150];
		sprintf(line,"'%s' using 1:%u title '%s' with linespoints lw 5 ps 2%s\n",
					 out_file.c_str(),i+2,param_names[i],tail);
		o << line;
	}
	o.close();

	system( ("gnuplot "+gnuplot_file).c_str() );
}

