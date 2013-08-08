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


class CFLSteadyAnalyser : public SequenceSteadyAnalyser<rheolef::Float>
{
public:

	CFLSteadyAnalyser( XMLConfigFile const& conf ):
		SequenceSteadyAnalyser<rheolef::Float>(
				conf.atoi("n_points_to_monitor"),
				conf.atof("normalized_deviation_limit") )
	{}
};

#endif /* CFLSTEADYANALYSER_H_ */
