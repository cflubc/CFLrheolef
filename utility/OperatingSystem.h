/*
 * OperatingSystem.h
 *
 *  Created on: 2013-05-23
 *      Author: ali
 */

#ifndef OPERATINGSYSTEM_H_
#define OPERATINGSYSTEM_H_


#include <string>

namespace OS {

void run_command( std::string const& command );
void mkdir( std::string const& directory );
void changedir( std::string const& directory );
std::string const& working_directory();

inline void
chdir_up()
{changedir("..");}

}

#endif /* OPERATINGSYSTEM_H_ */
