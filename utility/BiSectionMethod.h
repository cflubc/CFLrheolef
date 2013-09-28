/*
 * BiSectionMethod.h
 *
 *  Created on: Sep 24, 2013
 *      Author: ali
 */

#ifndef BISECTIONMETHOD_H_
#define BISECTIONMETHOD_H_

#include <cstddef>
#include <cmath>
#include <stdexcept>


struct BiSection
{
	template< typename T, typename Function >
	static T solve( T const& a, T const& b, int const niter, Function const& f )
	{
		check_sings(a,b,f);

		T l1 = a;
		T l2 = b;
		for(int i=0; i<niter; ++i)
		{
			T const mid = (l1+l2)/2.;
			if( f(l1)*f(mid)<0 )
				l2 = mid;
			else
				l1 = mid;
		}
		return (l1+l2)/2.;
	}

	template< typename T, typename Function >
	static T solve( T const& a, T const& b, T const& tol, Function const& f )
	{
		using namespace std;
		check_signs(a,b,f);

		T const diff = abs(b-a);
		if( diff<tol )
			return (a+b)/2.;
		int const niter = ceil( log2(diff/tol) );
		solve(a,b,niter,f);
	}

private:

	template< typename T, typename Function >
	static void check_sings( T const& a, T const& b, Function const& f ){
		if( 0<f(a)*f(b) )
			throw std::logic_error("Bisection method inputs shouldn't have the same sign");
	}
};



#endif /* BISECTIONMETHOD_H_ */
