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

#include "NormalStressBC_RHS.h"
#include "BodyForce.h"

#include "BubbleEncapsulationMesh.h"
#include "ChannelMesh.h"
#include "WavyChannelMesh.h"


struct voidMesh {
	voidMesh( XMLConfigFile const&, std::string const& ) {}
};

struct VoidRHS {
	VoidRHS( XMLConfigFile const&, rheolef::space const& ) {}
	void add_to_rhs( rheolef::field& ) const {}
};

typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> StokesFlow;


struct Problem_NewtonianCavity
{
	typedef StokesFlow Application;
	typedef FlowFields FieldsPool;
	typedef cavityBC BC;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "NewtonianCavity";
};

struct Problem_WavyChannelFouling
{
	typedef AugmentedLagrangianUnitFlow<StokesFlow,BodyForce> Application;
	typedef channelBC  BC;
	typedef FlowFields FieldsPool;
	typedef WavyChannelMesh Mesh;
	static constexpr cstr Name = "WavyChannelFouling";
};

struct Problem_AugLag_ChannelUnitFlow
{
	typedef AugmentedLagrangianUnitFlow<StokesFlow,BodyForce> Application;
	typedef channelBC  BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_ChannelUnitFlow";
};

struct Problem_AugLag_SteadyPoiseuille
{
	typedef StandardAugmentedLagrangian<StokesFlow,BodyForce> Application;
	typedef channelBC  BC;
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
	typedef bubble_BC  BC;
	typedef FlowFields FieldsPool;
	typedef BubbleEncapsulationMesh Mesh;
	static constexpr cstr Name = "AugLag_BubbleEncapsulation";
};



#endif /* PROBLEMS_H_ */
