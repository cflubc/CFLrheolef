/*
 * PrintArguments.h
 *
 *  Created on: 2013-05-23
 *      Author: ali
 */

#ifndef PRINTARGUMENTS_H_
#define PRINTARGUMENTS_H_

#include <utility>


template< typename Stream >
inline void
print_args( Stream& s )
{}

template< typename Stream, typename T, typename... Args >
inline void
print_args( Stream& s, T&& t, Args&&... args )
{
	s << std::forward<T>(t);
	print_args(s, std::forward<Args>(args)...);
}


template< typename Stream, typename... Args >
inline void
println_args( Stream& s, Args&&... args )
{
	print_args(s, std::forward<Args>(args)..., '\n');
}


template< typename Stream >
inline void
extract_args( Stream& s )
{}

template< typename Stream, typename T, typename... Args >
inline void
extract_args( Stream& s, T& t, Args&... args )
{
	s >> t;
	extract_args(s, args...);
}


#endif /* PRINTARGUMENTS_H_ */
