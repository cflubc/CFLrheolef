/*
 * OutputFormatting.h
 *
 *  Created on: 2013-04-29
 *      Author: ali
 */

#ifndef OUTPUTFORMATTING_H_
#define OUTPUTFORMATTING_H_

#include <stdexcept>
#include <initializer_list>

/**
 * print a table with given column widths easily.
 */
template< int Ncolumns, typename Stream >
class ColumnOutputFormatter
{
public:
	typedef const std::initializer_list<int>  width_list;

	ColumnOutputFormatter( Stream& ostream, width_list columns_width );
	void fill_horizontal( char c );

	template< typename... Args >
	void print( const Args... args ){
		static_assert( sizeof...(Args)==Ncolumns, "Wrong number of inputs provided" );
		recurse_print(args...);
	}


private:

	Stream& out;
	int width[Ncolumns];
	int total_width;

	template< typename T, typename... Args >
	void recurse_print( const T head, const Args... tail){
		out.width( width[Ncolumns-sizeof...(Args)-1] );
		out << head;
		recurse_print(tail...);
	}

	void recurse_print() const
	{out<<'\n';}
};



template< int Ncolumns, typename Stream  >
ColumnOutputFormatter<Ncolumns,Stream>::ColumnOutputFormatter(
			Stream& ostream,
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


template< int Ncolumns, typename Stream  >
void ColumnOutputFormatter<Ncolumns,Stream>::fill_horizontal( char c )
{
	const char old_fill = out.fill();
	out.fill(c);
	out.width(total_width);
	out << '\n';
	out.fill(old_fill);
}


template< typename Stream, typename... WidthList >
inline ColumnOutputFormatter<sizeof...(WidthList),Stream>
make_column_output( Stream& out, WidthList... list )
{
	return ColumnOutputFormatter<sizeof...(WidthList),Stream>(out,{list...});
}



#endif /* OUTPUTFORMATTING_H_ */
