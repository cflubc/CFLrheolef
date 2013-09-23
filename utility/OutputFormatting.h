/*
 * OutputFormatting.h
 *
 *  Created on: 2013-04-29
 *      Author: ali
 */

#ifndef OUTPUTFORMATTING_H_
#define OUTPUTFORMATTING_H_

#include <utility>
#include <numeric>
#include <stdexcept>

#include "PrintArguments.h"


/**
 * print a table with given column widths easily.
 */
template< int Ncolumns, typename Stream >
class ColumnOutputFormatter
{
public:

	template< typename... Nums >
	ColumnOutputFormatter( Stream& ostream, Nums&&... widths );

	void fill_horizontal( char c );

	template< typename... Args >
	void print( Args&&... args ){
		static_assert( sizeof...(Args)==Ncolumns, "Wrong number of inputs provided" );
		recurse_print( std::forward<Args>(args)...);
	}

	int table_width() const
	{return total_width;}

	void newline() const
	{out<<'\n';}

	template< typename... Args >
	void unformatted_print( Args&&... args ){
		print_args(out, std::forward<Args>(args)...);
	}

private:

	Stream& out;
	int const width[Ncolumns];
	int const total_width;

	template< typename T, typename... Args >
	void recurse_print( T&& head, Args&&... tail){
		out.width( width[Ncolumns-sizeof...(Args)-1] );
		out << head;
		recurse_print( std::forward<Args>(tail)... );
	}

	void recurse_print() const
	{newline();}
};


template< int Ncolumns, typename Stream >
template< typename... Nums >
ColumnOutputFormatter<Ncolumns,Stream>::ColumnOutputFormatter(
			Stream& ostream,
			Nums&&... widths ):
	out(ostream),
	width{widths...},
	total_width( std::accumulate(std::begin(width), std::end(width), 0) )
{}


template< int Ncolumns, typename Stream  >
void ColumnOutputFormatter<Ncolumns,Stream>::fill_horizontal( char c )
{

	const char old_fill = out.fill();
	out.fill(c);
	out.width(total_width);
	newline();
	out.fill(old_fill);
}


template< typename Stream, typename... WidthList >
inline ColumnOutputFormatter<sizeof...(WidthList),Stream>
make_column_output( Stream& out, WidthList&&... list )
{
	return ColumnOutputFormatter<sizeof...(WidthList),Stream>(out,std::forward<WidthList>(list)...);
}



#endif /* OUTPUTFORMATTING_H_ */
