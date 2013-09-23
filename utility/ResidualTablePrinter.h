/*
 * ResidualTablePrinter.h
 *
 *  Created on: 2013-06-02
 *      Author: ali
 */

#ifndef RESIDUALTABLEPRINTER_H_
#define RESIDUALTABLEPRINTER_H_

#include <cstddef>
#include <utility>
#include <string>
#include <stdexcept>

#include "RecuringAlarm.h"
#include "OutputFormatting.h"


template< int Ncolumns, typename Stream >
class ResidualTablePrinter
{
	ColumnOutputFormatter<Ncolumns,Stream> table;
	RecuringAlarm header_reprint;
	std::string const dashed_line;

public:

	template< typename... WidthList >
	ResidualTablePrinter(
							std::size_t const header_reprint_frequency,
							Stream& s,
							WidthList&& ... list ):
		table(s, std::forward<WidthList>(list)...),
		header_reprint(header_reprint_frequency),
		dashed_line(table.table_width(),'-')
	{}

	template< typename... Args >
	void print( Args&&... args )
	{table.print( std::forward<Args>(args)... );}

	template< typename... Args >
	void print_header_if_needed( Args&&... args )
	{
		if( header_reprint.alarm_ringing() ){
			table.newline();
			table.print( std::forward<Args>(args)... );
			table.unformatted_print(dashed_line);
			table.newline();
		}
	}
};


template< typename Stream, typename... WidthList >
inline ResidualTablePrinter<sizeof...(WidthList),Stream>
make_residual_table(
		std::size_t const header_reprint_frequency,
		Stream& out,
		WidthList&&... list )
{
	return ResidualTablePrinter<sizeof...(WidthList),Stream>
	        (header_reprint_frequency, out, std::forward<WidthList>(list)...);
}


#endif /* RESIDUALTABLEPRINTER_H_ */
