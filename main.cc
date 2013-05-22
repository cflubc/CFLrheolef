/*
 * stokes.cc
 *
 *  Created on: 2013-04-07
 *      Author: ali
 */

#include <cstdlib>
#include <string>
#include "rheolef.h"
#include "rheolef/diststream.h"

#include "ConfigXML.h"
#include "BlockSystem_abtb.h"
#include "BCs.h"
#include "IncompressibleStokesSolver.h"
#include "adaptationLoop.h"
#include "adaptationCriterions.h"
#include "StandardAugmentedLagrangian.h"
#include "DiffusionForms.h"
#include "AugmentedLagrangianUnitFlow.h"

// saving for restart: is a solver adapter? or new wrapper object?
/*
 * restart class: it runs the whole system first, if it is a case
 * of restart, reads the proper geometry and pass to adapt loop
 */

typedef FlowFields FieldsPool;
typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> StokesFlow;


struct Problem_ChannelUnitFlow
{
	typedef AugmentedLagrangianUnitFlow<StokesFlow,NormalStressBC_RHS> Application; //voidrhs
	typedef channel_fullBC  BC;
};

struct Problem_AugmentedLagrangian_SteadyPoiseuille
{
	typedef StandardAugmentedLagrangian<StokesFlow,NormalStressBC> Application;
	typedef channel_fullBC  BC;
};

struct Problem_AugmentedLagrangian_SteadyCavity
{
	typedef StandardAugmentedLagrangian<StokesFlow,VoidRHS> Application;
	typedef cavityBC BC;
};


//typedef AdaptationLoop<SApplication,DirichletBoundaryConditions>  Application;
typedef Problem_ChannelUnitFlow Problem;
typedef Problem::BC DirichletBoundaryConditions;
typedef Problem::Application  Application;

int main(int argc, char** argv )
{
	rheolef::environment env(argc,argv);
	rheolef::geo omega( argv[1] );
	const XMLConfigFile conf( argv[2] );

	const DirichletBoundaryConditions BC( conf.child("BC") );
	FieldsPool fields(conf,omega,BC);

	Application app(conf,fields,BC);
	app.run();

	return EXIT_SUCCESS;
}


template< typename Application >
struct NonAdaptiveApplication
{
	template< typename Fields, typename BoundaryCondition >
	void setup_and_run( const XMLConfigFile& conf, Fields& fields, BoundaryCondition& BC ){
		Application app(conf,fields);
		app.run();
	}
};
