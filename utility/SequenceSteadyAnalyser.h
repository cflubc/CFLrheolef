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

	SequenceSteadyAnalyser( int const n, T const& limit ):
		nsample(n),
		steady_state_limit(limit)
	{
		if(n<3)
			throw std::logic_error("Sequence analysis needs at leat 3 points");
	}

	template< typename Iterator >
	bool sequence_steady_state_reached( Iterator first, Iterator last );

	template< typename Container >
	bool sequence_steady_state_reached( Container const& C )
	{return sequence_steady_state_reached( C.begin(), C.end() );}

private:

	T const nsample;
	T const steady_state_limit;
};


template< typename T >
template< typename Iterator >
bool SequenceSteadyAnalyser<T>::sequence_steady_state_reached( Iterator first, Iterator last )
{
	if( last-first<nsample )
		return false;

	first = last-nsample;
	T const ave = std::accumulate(first,last,0.)/nsample;
	T normalized_deviation = 0.;
	for( Iterator i=first; i!=last; ++i )
		normalized_deviation += fabs(*i-ave);
	normalized_deviation /= nsample*ave;

	return (normalized_deviation < steady_state_limit);
}

#endif /* SEQUENCESTEADYANALYSER_H_ */
