/*
 * TensorFieldIterator.h
 *
 *  Created on: 2013-05-01
 *      Author: ali
 */

#ifndef TENSORFIELDITERATOR_H_
#define TENSORFIELDITERATOR_H_

#include "rheolef/compiler.h"


template< typename Iterator >
class TensorFieldIterator
{
	typedef rheolef::Float Float;

public:
	enum { Ncomp=3 };
	typedef Float tensor[Ncomp];

	/// for norm of symmetric tensor: @f$ T^2_{xx} + 2T^2_{xy} + T^2_{yy} @f$
	constexpr static tensor coef_for_norm_calc{1., 2., 1.};

	template< typename Field >
	TensorFieldIterator( Field& f ):
		it{ f[0].begin_dof(),
		    f[1].begin_dof(),
		    f[2].begin_dof() },
		end0( f[0].end_dof() )
	{}

	TensorFieldIterator() {}

	typename Iterator::reference operator() ( int i ) const
	{ return *(it[i]); }

	void operator++(){
		++(it[0]);
		++(it[1]);
		++(it[2]);
	}

	bool end_reached() const
	{ return it[0]==end0; }


private:
	Iterator it[Ncomp];
	Iterator end0;
};


template< typename Iterator >
constexpr typename TensorFieldIterator<Iterator>::tensor
TensorFieldIterator<Iterator>::coef_for_norm_calc;


#endif /* TENSORFIELDITERATOR_H_ */
