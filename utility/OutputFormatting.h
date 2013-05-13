/*
 * OutputFormatting.h
 *
 *  Created on: 2013-04-29
 *      Author: ali
 */

#ifndef OUTPUTFORMATTING_H_
#define OUTPUTFORMATTING_H_

#include <cmath>
#include <string>
#include <ostream>
#include <stdexcept>
#include <initializer_list>

/**
 * print a table with given column widths easily.
 */
template< int Ncolumns >
class ColumnOutputFormatter
{
	typedef const std::initializer_list<int>  width_list;

public:
	ColumnOutputFormatter( std::ostream& ostream, width_list columns_width );
	void fill_horizontal( char c, double extend_coef=1. );

	template< typename... Args >
	void print( const Args... args ){
		static_assert( sizeof...(Args)==Ncolumns, "Wrong number of inputs provided" );
		recurse_print(args...);
	}


private:
	std::ostream& out;
	int width[Ncolumns];
	int total_width;

	template< typename T, typename... Args >
	void recurse_print( const T head, const Args... tail){
		out.width( width[Ncolumns-sizeof...(Args)-1] );
		out << head;
		recurse_print(tail...);
	}

	void recurse_print()
	{out<<'\n';}
};



template< int Ncolumns >
ColumnOutputFormatter<Ncolumns>::ColumnOutputFormatter(
			std::ostream& ostream,
			width_list columns_width ):
	out(ostream)
{
	total_width = 0;
	int i=0;
	for(auto x : columns_width){
		width[i] = x;
		total_width += x;
		++i;
	}
}


template< int Ncolumns >
void ColumnOutputFormatter<Ncolumns>::fill_horizontal( char c, double extend_coef )
{
	const char old_fill = out.fill();
	out.fill(c);
	out.width( std::ceil(extend_coef*total_width) );
	out << '\n';
	out.fill(old_fill);
}


template< typename... WidthList >
inline ColumnOutputFormatter<sizeof...(WidthList)>
make_column_output( std::ostream& out, WidthList... list )
{
	return ColumnOutputFormatter<sizeof...(WidthList)>(out,{list...});
}



#endif /* OUTPUTFORMATTING_H_ */
