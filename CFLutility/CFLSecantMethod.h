/*
 * CFLSecantMethod.h
 *
 *  Created on: 2013-05-21
 *      Author: ali
 */

#ifndef CFLSECANTMETHOD_H_
#define CFLSECANTMETHOD_H_

#include <cstddef>
#include <sstream>
#include "rheolef.h"

#include "ConfigXML.h"
#include "SecantMethod.h"



class CFLSecantMethod : public SecantMethod<rheolef::Float>
{
	typedef rheolef::Float Float;

public:

	CFLSecantMethod( XMLConfigFile const& conf, Float const& target ):
		SecantMethod<Float>(
				conf("max_iter",std::size_t()),
				conf("tolerance",Float()),
				target,
				x1(init_point_str(conf)),
				f1(init_point_str(conf)),
				conf("next_input",Float())
			  	  	  	  )
	{}

private:
	static cstr init_point_str( XMLConfigFile const& cf )
	{return cf("initial_point");}

	static Float x1( cstr const str ){
		std::istringstream s(str);
		Float x1;
		s >> x1;
		return x1;
	}

	static Float f1( cstr const str ){
		std::istringstream s(str);
		Float f1;
		s >> f1;
		s >> f1;
		return f1;
	}
};


#endif /* CFLSECANTMETHOD_H_ */
