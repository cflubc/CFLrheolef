/*
 * Problems.cc
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */

#include "Problems.h"

#define DEFINE_IN_NAMESPACE(problem) \
		constexpr cstr problem::Name;


DEFINE_IN_NAMESPACE(Problem_NewtonianCavity)
DEFINE_IN_NAMESPACE(Problem_WavyChannelFouling)
DEFINE_IN_NAMESPACE(Problem_AugLag_ChannelUnitFlow)
DEFINE_IN_NAMESPACE(Problem_AugLag_SteadyPoiseuille)
DEFINE_IN_NAMESPACE(Problem_AugLag_SteadyCavity)
DEFINE_IN_NAMESPACE(Problem_AugLag_BubbleEncapsulation)
