/*
 * Problems.h
 *
 *  Created on: 2013-05-26
 *      Author: ali
 */

#ifndef PROBLEMS_H_
#define PROBLEMS_H_

#include "CFL.h"
#include "BCs.h"
#include "DiffusionForms.h"
#include "BlockSystem_abtb.h"
#include "IncompressibleStokesSolver.h"
#include "StandardAugmentedLagrangian.h"
#include "AugmentedLagrangianUnitFlow.h"

#include "BubbleEncapsulationMesh.h"
#include "ChannelMesh.h"


struct voidMesh
{
	voidMesh( XMLConfigFile const&, std::string const& )
	{}
};


typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> StokesFlow;


struct Problem_AugLag_ChannelUnitFlow
{
	typedef AugmentedLagrangianUnitFlow<StokesFlow,NormalStressBC_RHS> Application;
	typedef channel_fullBC  BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_ChannelUnitFlow";
};

struct Problem_AugLag_SteadyPoiseuille
{
	typedef StandardAugmentedLagrangian<StokesFlow,NormalStressBC_RHS> Application;
	typedef channel_fullBC  BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_SteadyPoiseuille";
};

struct Problem_AugLag_SteadyCavity
{
	typedef StandardAugmentedLagrangian<StokesFlow,VoidRHS> Application;
	typedef cavityBC BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_SteadyCavity";
};


struct Problem_AugLag_BubbleEncapsulation
{
	typedef StandardAugmentedLagrangian<StokesFlow,NormalStressBC_RHS> Application;
	typedef channel_fullBC  BC;
	typedef FlowFields FieldsPool;
	typedef BubbleEncapsulationMesh Mesh;
	static constexpr cstr Name = "AugLag_BubbleEncapsulation";
};



#endif /* PROBLEMS_H_ */
