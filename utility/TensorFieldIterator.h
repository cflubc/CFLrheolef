/*
 * TensorFieldIterator.h
 *
 *  Created on: 2013-05-01
 *      Author: ali
 */

#ifndef TENSORFIELDITERATOR_H_
#define TENSORFIELDITERATOR_H_

#include "rheolef.h"


template< typename Iterator >
class TensorFieldIterator
{
	typedef rheolef::field field;
	typedef rheolef::Float Float;

public:
	enum { Ncomp=3 };
	typedef Float tensor[Ncomp];

	/// for norm of symmetric tensor: @f$ T^2_{xx} + 2T^2_{xy} + T^2_{yy} @f$
	constexpr static tensor coef_for_norm_calc{1., 2., 1.};

	TensorFieldIterator( field& f ):
		it{ f[0].begin_dof(),
		    f[1].begin_dof(),
		    f[2].begin_dof() },
		end0( f[0].end_dof() )
	{}

	TensorFieldIterator() {}

	typename Iterator::reference operator() ( int i ) const
	{ return *(it[i]); }

	void operator++(){
//		for(int i=0; i<Ncomp; ++i)
//			++(it[i]);
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



//template< typename Iterator >
//class TensorFieldIterator
//{
//	typedef rheolef::field field;
//
//public:
//	enum { Ncomp=4 };
//	typedef rheolef::Float tensor[Ncomp];
//
//	TensorFieldIterator( field& f ):
//		it{ f(0,0).begin_dof(),
//		    f(1,0).begin_dof(),
//		    f(0,1).begin_dof(),
//		    f(1,1).begin_dof() },
//		end00(f(0,0).end_dof())
//	{}
//
//	TensorFieldIterator() {}
//
//	typename Iterator::reference operator() ( int i ) const
//	{ return *(it[i]); }
//
//	void operator++(){
//		for(int i=0; i<Ncomp; ++i)
//			++(it[i]);
//	}
//
//	bool end_reached()
//	{ return it[0]==end00; }
//
//	void reinit( field& f ){
//		it[0] = f(0,0).begin_dof();
//		it[1] = f(1,0).begin_dof();
//		it[2] = f(0,1).begin_dof();
//		it[3] = f(1,1).begin_dof();
//	}
//
//private:
//	Iterator it[Ncomp];
//	Iterator end00;
//};
