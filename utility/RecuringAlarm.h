/*
 * RecuringAlarm.h
 *
 *  Created on: 2013-06-03
 *      Author: ali
 */

#ifndef RECURINGALARM_H_
#define RECURINGALARM_H_

#include <cstddef>

class RecuringAlarm
{
	typedef std::size_t size_t;

	size_t const n_recuring;
	size_t counter;

public:

	RecuringAlarm( size_t const n, size_t const initval=0 ):
		n_recuring(n),
		counter(initval)
	{}

	bool alarm_ringing()
	{ return (counter++%n_recuring)==0; }

	void reset()
	{ counter=0; }
};



#endif /* RECURINGALARM_H_ */
