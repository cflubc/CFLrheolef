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
#include "SteadyStokesAugmentedIteration.h"
#include "AugmentedLagrangianUnitFlow.h"
#include "IncompressibleNavierStokes.h"
#include "SteadyNavierStokesAugmentedLagrangian.h"

#include "NormalStressBC_RHS.h"
#include "BodyForce.h"

#include "BubbleEncapsulationMesh.h"
#include "ChannelMesh.h"
#include "WavyChannelMesh.h"
#include "FlowOnsetDetection.h"
#include "WavyFracture.h"


struct VoidRHS {
	VoidRHS( XMLConfigFile const&, rheolef::space const& ) {}
	void add_to_rhs( rheolef::field& ) const {}
};

typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> StokesFlow;
typedef AugmentedLagrangian_basic<AugmentedLagrangian_param::Bingham_unique,AugmentedLagrangian_param::alpha_unique> ALbasic_unique_params;
typedef AugmentedLagrangian_basic<AugmentedLagrangian_param::Bingham_multiRegion,AugmentedLagrangian_param::alpha_unique> ALbasic_multiRegion_Bn;
typedef AugmentedLagrangian_basic<AugmentedLagrangian_param::Bingham_multiRegion,AugmentedLagrangian_param::alpha_multiRegion> ALbasic_multiRegion;

struct Problem_NewtonianCavity
{
	typedef StokesFlow Application;
	typedef FlowFields FieldsPool;
	typedef cavityBC BC;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "NewtonianCavity";
};

struct Problem_NavierStokesCavity
{
//	typedef SteadyNavierStokesAugmentedLagrangian<ALbasic_unique_params> Application;
//	typedef IncompressibleNavierStokes Application;
	typedef SteadyNavierStokesAugmentedLagrangian<ALbasic_unique_params,VoidRHS> Application;
	typedef FlowFields FieldsPool;
	typedef cavityBC BC;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "NavierStokesCavity";
};

struct Problem_WavyFouling
{
//	typedef AugmentedLagrangianUnitFlow<ALbasic_unique_params,BodyForce> Application;
	typedef SteadyNavierStokesAugmentedLagrangian<ALbasic_unique_params,BodyForce> Application;
	typedef channelBC BC;
	typedef FlowFields FieldsPool;
	typedef WavyChannelMesh Mesh;
	static constexpr cstr Name = "WavyFouling";
};

struct Problem_AugLag_ChannelUnitFlow
{
	typedef AugmentedLagrangianUnitFlow<ALbasic_unique_params,BodyForce> Application;
	typedef channelBC  BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_ChannelUnitFlow";
};


struct Problem_AugLag_SteadyPoiseuille
{
	typedef SteadyStokesAugmentedLagrangian<ALbasic_unique_params,BodyForce> Application;
	typedef channelBC  BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_SteadyPoiseuille";
};

struct Problem_AugLag_SteadyCavity
{
	typedef SteadyNavierStokesAugmentedLagrangian<ALbasic_unique_params,VoidRHS> Application;
	typedef cavityBC BC;
	typedef FlowFields FieldsPool;
	typedef ChannelMesh Mesh;
	static constexpr cstr Name = "AugLag_SteadyCavity";
};

struct Problem_BubbleEncapsulation
{
	typedef AugmentedLagrangianUnitFlow<ALbasic_unique_params,NormalStressBC_RHS> Application;
	typedef channelBC  BC;
	typedef FlowFields FieldsPool;
	typedef BubbleEncapsulationMesh Mesh;
	static constexpr cstr Name = "AugLag_BubbleEncapsulation";
};

struct Problem_DropletEncapsulation
{
	typedef AugmentedLagrangianUnitFlow<ALbasic_multiRegion_Bn,BodyForce> Application;
	typedef channelBC BC;
	typedef FlowFields FieldsPool;
	typedef BubbleEncapsulationMesh Mesh;
	static constexpr cstr Name = "AugLag_DropletEncapsulation";
};

struct Problem_MacroBubbleFlowOnset
{
	typedef FlowOnsetDetection<ALbasic_unique_params,NormalStressBC_RHS> Application;
	typedef bubble_BC  BC;
	typedef FlowFields FieldsPool;
	typedef BubbleEncapsulationMesh Mesh;
	static constexpr cstr Name = "MacroBubbleFlowOnset";
};

#endif /* PROBLEMS_H_ */
