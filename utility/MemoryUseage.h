/**
 * MemoryUseage.h
 *
 *  Created on: 2013-05-22
 *      Author: ali
 *      
 *  The code here is inspired from Pelicans library
 */

#ifndef MEMORYUSEAGE_H_
#define MEMORYUSEAGE_H_

#include <cstdlib>
#include "PrintArguments.h"

double
retrieve_size_in_status( std::string const& keyword );


inline double
used_memory()
{return retrieve_size_in_status("VmRSS:");}

inline double
peak_used_memory()
{return retrieve_size_in_status("VmHWM:");}


template< typename stream >
void print_memory_useage( stream& s)
{
	print_args(s,"Memory useage:\n",
		         "  Average: ", used_memory(),      " Mb\n",
		         "     Peak: ", peak_used_memory(), " Mb\n"
			    );
}

#endif /* MEMORYUSEAGE_H_ */
