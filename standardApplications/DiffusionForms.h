/*
 * DiffusionForms.h
 *
 *  Created on: 2013-05-06
 *      Author: ali
 */

#ifndef DIFFUSIONFORMS_H_
#define DIFFUSIONFORMS_H_

#include "rheolef.h"
#include "ConfigXML.h"


struct ConstantDiffusionForm
{
	ConstantDiffusionForm( const XMLConfigFile& conf,
			 	 	 	   const rheolef::space& Xh ):
		D(Xh,Xh,"2D_D")
	{}

	ConstantDiffusionForm( const XMLConfigFile& conf,
			               const rheolef::space& Xh,
			               const rheolef::Float& r ):
		ConstantDiffusionForm(conf,Xh)
	{D*=r;}

	ConstantDiffusionForm()
	{}

	rheolef::form get_form() const
	{return D;}

	rheolef::form D;
};


#endif /* DIFFUSIONFORMS_H_ */

