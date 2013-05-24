/*
 * OperatingSystem.cc
 *
 *  Created on: 2013-05-23
 *      Author: ali
 */


#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <cstdlib>
#include "OperatingSystem.h"

namespace OS {

void mkdir( std::string const& directory )
{
   DIR* dir = opendir( directory.c_str() );
   if( dir!=0 )
	  closedir(dir);
   else
   {
	  std::string line;
	  line = "mkdir "+directory;
	  int sys_res = system( line.c_str() ) ;
   }
}


void changedir( std::string const& directory )
{
	chdir( directory.c_str() );
}


std::string const& working_directory()
{
	char tmp[1024];
	getcwd(tmp,sizeof(tmp));
	std::string cwd(tmp);
	return cwd;
}


}
