/**
 * CFLSteadyAnalyser.h
 *
 *  Created on: 2013-05-29
 *      Author: ali
 *
 * This is only a wrapper class for SequenceSteadyAnalyser class.
 * In fact it merely adds a simple constructor to initialize from
 * a XML config file.
 */

#ifndef CFLSTEADYANALYSER_H_
#define CFLSTEADYANALYSER_H_

#include "rheolef/compiler.h"

#include "ConfigXML.h"
#include "SequenceSteadyAnalyser.h"


class CFLSteadyAnalyser
{
	SequenceSteadyAnalyser<rheolef::Float> seq;

public:

	CFLSteadyAnalyser( XMLConfigFile const& conf ):
		seq( conf.atoi("n_points_to_monitor"), conf.atof("normalized_deviation_limit") )
	{}

	template< typename Container >
	bool sequence_steady_state_reached( Container const& C )
	{return seq.sequence_steady_state_reached(C);}
};

#endif /* CFLSTEADYANALYSER_H_ */
