/*
 * TimeGauge.h
 *
 *  Created on: 2013-05-22
 *      Author: ali
 */

#ifndef TIMEGAUGE_H_
#define TIMEGAUGE_H_

#include <chrono>
#include <ratio>

class TimeGauge
{
	typedef std::chrono::steady_clock clock;
	typedef clock::time_point time_point;
	typedef std::chrono::duration<double,std::ratio<3600>> hours;

public:
	void start()
	{tbegin = clock::now();}

	void stop()
	{tend = clock::now();}

	hours::rep get_time_passed() const {
		 return std::chrono::duration_cast<hours>(tend-tbegin).count();
	}

private:
	time_point tbegin;
	time_point tend;
};


#endif /* TIMEGAUGE_H_ */
