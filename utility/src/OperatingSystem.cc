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
#include <cstdio>
#include "OperatingSystem.h"

namespace OS {


void run_command( std::string const& command )
{
	printf( "[Execute] %s\n", command.c_str() );
	int sys_res = system( command.c_str() );
}


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
