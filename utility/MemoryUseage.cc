/*
 * MemoryUseage.cc
 *
 *  Created on: 2013-05-22
 *      Author: ali
 */



#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

// Linux specific headers
#include <sys/types.h>
#include <unistd.h>

#include "MemoryUseage.h"


double
retrieve_size_in_status( std::string const& keyword )
{
   std::ostringstream os;
   os << "/proc/" << getpid() << "/status";
   std::ifstream in( os.str().c_str() ) ;
   if(!in){
      printf( "Unable to open %s\n",os.str().c_str() );
      return 0;
   }


  double result = 0;
  while( !in.eof() )
  {
	 std::string word;
	 in >> word;
	 if( word==keyword )
	 {
		in >> result ;
		in >> word ;
		//conver to Mb
		if( word=="kB" )
			result /= 1024;
		else
			printf("Memory unit is %s instead of kB\n",word.c_str());
		break;
	 }
  }
  in.close();
  return result;
}


