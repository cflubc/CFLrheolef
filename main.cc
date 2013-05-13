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

// saving for restart: is a solver adapter? or new wrapper object?
/*
 * restart class: it runs the whole system first, if it is a case
 * of restart, reads the proper geometry and pass to adapt loop
 */

typedef FlowFields FieldsPool;

typedef IncompLinearDiffusionStokesSolver<BlockSystem_abtb> StokesFlow;
typedef StandardAugmentedLagrangian<StokesFlow> SApplication;

typedef SApplication Application;
typedef cavityBC  DirichletBoundaryConditions;
//typedef AdaptationLoop<SApplication,DirichletBoundaryConditions>  Application;


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
