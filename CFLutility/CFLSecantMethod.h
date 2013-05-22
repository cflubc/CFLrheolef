/*
 * CFLSecantMethod.h
 *
 *  Created on: 2013-05-21
 *      Author: ali
 */

#ifndef CFLSECANTMETHOD_H_
#define CFLSECANTMETHOD_H_

#include <sstream>
#include "rheolef.h"

#include "ConfigXML.h"
#include "SecantMethod.h"



class CFLSecantMethod
{
	typedef rheolef::Float Float;

public:
	CFLSecantMethod( XMLConfigFile const& conf ):
		secant( conf.atoi("max_iter"),
				conf.atof("tolerance"),
				conf.atof("target"),
				x1(init_point_str(conf)),
				f1(init_point_str(conf)),
				conf.atof("next_input")
			  )
	{}

	void reset()
	{secant.reset();}

	size_t n_iterations_done() const
	{return secant.n_iterations_done();}

	bool not_converged_and_have_iterations_left() const
	{return secant.not_converged_and_have_iterations_left();}

	Float get_input() const
	{return secant.get_input();}

	Float predict_new_input( Float const& f2 )
	{return secant.predict_new_input(f2);}

private:

	struct Float_abs {
		static Float abs( Float const& x )
		{return rheolef::abs(x);}
	};
	SecantMethod<Float,Float_abs> secant;

	static cstr init_point_str( XMLConfigFile const& cf )
	{return cf("initial_point");}

	static Float x1( cstr str ){
		std::istringstream s(str);
		Float x1;
		s >> x1;
		return x1;
	}

	static Float f1( cstr str ){
		std::istringstream s(str);
		Float f1;
		s >> f1;
		s >> f1;
		return f1;
	}
};


#endif /* CFLSECANTMETHOD_H_ */
