/*
 * MathUtility.h
 *
 *  Created on: 2013-06-27
 *      Author: ali
 */

#ifndef MATHUTILITY_H_
#define MATHUTILITY_H_


/** Options for include/exclude in RangeChecker class */
struct interval_constants
{
	constexpr static bool
		inclusive=true,   exclusive=false,
		include_beg=true, exclude_beg=false,
		include_end=true, exclude_end=false;
};

/**
 * Checks whether a number is in range beg,end. User instructs wether wants
 * to include/exclude end points in range.
 *
 * Examples:
 *    RangeChecker<double,range_constants::inclusive> x(0,1); //(0,1)
 *    RangeChecker<double,range_constants::include_beg,range_constants::exclude_end> x(0,1); //[0,1)
 */
template< typename T, bool including_beg, bool including_end = including_beg >
class Interval
{
public:

	Interval( T const& endpoint1, T const& endpoint2 ):
		was_initially_given_reversed( endpoint2<endpoint1 ),
		beg( was_initially_given_reversed ? endpoint2:endpoint1 ),
		end( was_initially_given_reversed ? endpoint1:endpoint2 )
	{}

	bool includes( T const& t ) const {
		return ( including_beg ? beg<=t:beg<t ) &&
			   ( including_end ? t<=end:t<end );
	}

	constexpr bool includes_begin() const
	{return including_beg;}

	constexpr bool includes_end() const
	{return including_end;}

	bool const was_initially_given_reversed;
	T const beg;
	T const end;
};

/** Typedefs for convinience */
template< typename T >
using ExclusiveInterval = Interval<T,interval_constants::exclusive>;

template< typename T >
using InclusiveInterval = Interval<T,interval_constants::inclusive>;

#endif /* MATHUTILITY_H_ */
