/*
 * PrintArguments.h
 *
 *  Created on: 2013-05-23
 *      Author: ali
 */

#ifndef PRINTARGUMENTS_H_
#define PRINTARGUMENTS_H_


template< typename Stream >
inline void
print_args( Stream& s )
{}

template< typename Stream, typename T, typename... Args >
inline void
print_args( Stream& s, T const& t, Args const&... args )
{
	s << t;
	print_args(s,args...);
}


template< typename Stream, typename... Args >
inline void
println_args( Stream& s, Args const&... args )
{
	print_args(s,args...,'\n');
}


#endif /* PRINTARGUMENTS_H_ */
