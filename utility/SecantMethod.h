/*
 * SecantMethod.h
 *
 *  Created on: 2013-05-20
 *      Author: ali
 */

#ifndef SECANTMETHOD_H_
#define SECANTMETHOD_H_

#include <cstdlib>
#include <cassert>
#include <cmath>


/**
 * Secant method for a function f with input parameter x.
 * @see en.wikipedia.org/wiki/Secant_method
 *
 * We have data of points 1,2 and predict the required x to reach
 * the target value using a linear approximation of function.
 */
template< typename T = double >
class SecantMethod
{
	T const target;
	T const tolerance;
	T f1;
	T x2;
	T dx;
	size_t max_iter;
	size_t iter;
	bool converged;

public:
	SecantMethod( size_t nmax, T const& tol, T const& _target, T const& x1, T const& _f1, T const& _x2 ):
		target(_target),
		tolerance(tol),
		f1(_f1),
		x2(_x2),
		dx(_x2-x1),
		max_iter(nmax),
		iter(0),
		converged(false)
	{}

	void reset(){
		iter = 0;
		converged = false;
	}

	size_t n_iterations_done() const
	{return iter;}

	bool not_converged_and_have_iterations_left() const
	{return (iter<max_iter) && !converged;}

	T get_input() const
	{return x2;}

	T predict_new_input( T const& f2 );
};


template< typename T >
T SecantMethod<T>::predict_new_input( T const& f2 )
{
	assert(f1!=f2);

	T const dtarget = target-f2;
	converged = fabs(dtarget)<tolerance;
	if( !converged ){
		++iter;
		dx = dx/(f2-f1)*dtarget;
		// update to new x
		x2 += dx;
		f1 = f2;
	}
	return x2;
}


#endif /* SECANTMETHOD_H_ */
