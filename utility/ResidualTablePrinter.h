/*
 * ResidualTablePrinter.h
 *
 *  Created on: 2013-06-02
 *      Author: ali
 */

#ifndef RESIDUALTABLEPRINTER_H_
#define RESIDUALTABLEPRINTER_H_

#include <cstddef>
#include <initializer_list>
#include <stdexcept>

#include "RecuringAlarm.h"
#include "OutputFormatting.h"


template< int Ncolumns, typename Stream >
class ResidualTablePrinter
{
	RecuringAlarm header_reprint;
	ColumnOutputFormatter<Ncolumns,Stream> table;

public:

	template< typename... WidthList >
	ResidualTablePrinter(
							std::size_t const header_reprint_frequency,
							Stream& s,
							WidthList... list ):
		header_reprint(header_reprint_frequency),
		table(s,{list...})
	{}

	template< typename... Args >
	void print( const Args... args )
	{table.print(args...);}

	template< typename... Args >
	void print_header_if_needed( const Args... args )
	{
		if( header_reprint.alarm_ringing() ){
			table.print(args...);
			table.fill_horizontal('-');
		}
	}
};


template< typename Stream, typename... WidthList >
inline ResidualTablePrinter<sizeof...(WidthList),Stream>
make_residual_table(
		std::size_t const header_reprint_frequency,
		Stream& out,
		WidthList... list )
{
	return ResidualTablePrinter<sizeof...(WidthList),Stream>(header_reprint_frequency,out,list...);
}


#endif /* RESIDUALTABLEPRINTER_H_ */
