/*
 * stokes.cc
 *
 *  Created on: 2013-04-07
 *      Author: ali
 */

#include <cstdlib>
#include <string>
#include <iomanip>
#include "rheolef.h"
#include "rheolef/diststream.h"

#include "CFL.h"
#include "ConfigXML.h"
#include "MemoryUseage.h"
#include "BlockSystem_abtb.h"
#include "BCs.h"
#include "IncompressibleStokesSolver.h"
#include "adaptationLoop.h"
#include "adaptationCriterions.h"
#include "StandardAugmentedLagrangian.h"
#include "DiffusionForms.h"
#include "AugmentedLagrangianUnitFlow.h"
#include "TimeGauge.h"
#include "PrintArguments.h"
#include "OperatingSystem.h"

// saving for restart: is a solver adapter? or new wrapper object?
/*
 * restart class: it runs the whole system first, if it is a case
 * of restart, reads the proper geometry and pass to adapt loop
 */

typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> StokesFlow;


struct Problem_ChannelUnitFlow
{
	typedef AugmentedLagrangianUnitFlow<StokesFlow,NormalStressBC_RHS> Application;
	typedef channel_fullBC  BC;
	typedef FlowFields FieldsPool;
};

struct Problem_AugmentedLagrangian_SteadyPoiseuille
{
	typedef StandardAugmentedLagrangian<StokesFlow,NormalStressBC_RHS> Application;
	typedef channel_fullBC  BC;
	typedef FlowFields FieldsPool;
};

struct Problem_AugmentedLagrangian_SteadyCavity
{
	typedef StandardAugmentedLagrangian<StokesFlow,VoidRHS> Application;
	typedef cavityBC BC;
	typedef FlowFields FieldsPool;
};


typedef Problem_AugmentedLagrangian_SteadyPoiseuille Problem;
typedef Problem::BC DirichletBoundaryConditions;
typedef Problem::FieldsPool FieldsPool;
//typedef Problem::Application  Application;
typedef AdaptationLoop<Problem::Application,Problem::FieldsPool,Problem::BC>  Application;

int main(int argc, char** argv )
{
	rheolef::environment env(argc,argv);
	rheolef::geo omega( argv[1] );
	XMLConfigFile const conf( argv[2] );

	TimeGauge timer;
	timer.start();

	XMLConfigFile const FE( conf.child(FieldsPool_Module) );
	DirichletBoundaryConditions BC(FE);
	FieldsPool fields(FE,omega,BC);

	CFL_mkresult_folder_and_cd_to_it(0);
	Application app(conf,fields,BC);
	app.run();

	timer.stop();

	rheolef::dout << std::setprecision(3);
	rheolef::dout << "\n--------------------------------\n";
	print_args(rheolef::dout,"Running time: ",timer.get_time_passed()," hours\n");
	print_memory_useage( rheolef::dout );
	rheolef::dout <<   "--------------------------------\n";

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
