/*
 * SequenceSteadyAnalyser.h
 *
 *  Created on: 2013-05-28
 *      Author: ali
 */

#ifndef SEQUENCESTEADYANALYSER_H_
#define SEQUENCESTEADYANALYSER_H_

#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>


template< typename T >
class SequenceSteadyAnalyser
{
public:

	SequenceSteadyAnalyser( int const n ):
		nsample(n),
		diff(n),
		first_time_sampling(true)
	{
		if(n<3)
			throw std::logic_error("Sequence analysis needs at leat 3 points");
	}

	template< typename Iterator >
	bool sequence_steady_state_reached( Iterator first, Iterator last )
	{
		if( last-first<nsample )
			return false;

		first = last-nsample;
		if( first_time_sampling ){
			auto res = std::minmax_element(first,last);
			sequence_min = *res.first;
			sequence_max = *res.second;
			first_time_sampling = false;
		}
		else {
			if( sequence_max<*last )
				sequence_max = *last;
			else
			if( *last<sequence_min )
				sequence_min = *last;
		}

		std::adjacent_difference( first, last, diff.begin() );
		sequence_absolute_diff_sum = 0.;
		for( auto i=diff.begin()+1; i!=diff.end(); ++i )
			sequence_absolute_diff_sum += fabs(*i);

		T const average_diff = sequence_absolute_diff_sum/(nsample-1);
		return (sequence_max-sequence_min)<steady_acceptance_coef*average_diff;
	}


private:

	int const nsample;
	T sequence_min;
	T sequence_max;
	T sequence_absolute_diff_sum;
	T steady_acceptance_coef;
	std::vector<T> diff;

	bool first_time_sampling;
};


#endif /* SEQUENCESTEADYANALYSER_H_ */
