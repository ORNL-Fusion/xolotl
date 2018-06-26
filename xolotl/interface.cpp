// Includes
#include "interface.h"
#include <iostream>
#include <Options.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <xolotlPerf.h>
//#include <INetworkLoader.h>
//#include <IReactionNetwork.h>
#include <SolverHandlerFactory.h>
//#include <ISolverHandler.h>
#include <IReactionHandlerFactory.h>

std::shared_ptr<xolotlSolver::PetscSolver> XolotlInterface::initializeXolotl() {
	// Initialize MPI with empty line
	char* mpiArgs[1];
	mpiArgs[0] = NULL;
	char** mpiFargv = mpiArgs;
	int mpiArgc = 0;
	MPI_Init(&mpiArgc, &mpiFargv);

	// Create the options
	xolotlCore::Options opts;

	// Create the command line
	string pathToFile("param.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();
	char* args[2];
	args[0] = const_cast<char*>(fname);
	args[1] = NULL;
	char** fargv = args;

	// Read the parameter file to set the options
	opts.readParams(fargv);

	// Set up our performance data infrastructure.
	xolotlPerf::initialize(opts.getPerfHandlerType());
	auto handlerRegistry = xolotlPerf::getHandlerRegistry();

	// Get the MPI rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Create the material factory
	auto materialFactory =
			xolotlFactory::IMaterialFactory::createMaterialFactory(
					opts.getMaterial(), opts.getDimensionNumber());
	// Initialize it with the options
	materialFactory->initializeMaterial(opts);

	// Initialize the temperature
	xolotlFactory::initializeTempHandler(opts);
	// Get the temperature handler
	auto tempHandler = xolotlFactory::getTemperatureHandler();

	// Initialize the visualization
	xolotlFactory::initializeVizHandler(opts.useVizStandardHandlers());

	// Create the network handler factory
	auto networkFactory =
				xolotlFactory::IReactionHandlerFactory::createNetworkFactory(
						opts.getMaterial());
	// Build a reaction network
	networkFactory->initializeReactionNetwork(opts, handlerRegistry);
	auto& network = networkFactory->getNetworkHandler();

	// Initialize and get the solver handler
	xolotlFactory::initializeDimension(opts, network);
	auto& solvHandler = xolotlFactory::getSolverHandler();
	// Initialize the solver handler
	solvHandler.initializeHandlers(materialFactory, tempHandler, opts);

	// Setup the solver
	std::shared_ptr<xolotlSolver::PetscSolver> solver(
			new xolotlSolver::PetscSolver(solvHandler, handlerRegistry));
	solver->setCommandLineOptions(opts.getPetscArgc(),
			opts.getPetscArgv());
	solver->initialize();

	return solver;
}

void XolotlInterface::solveXolotl(std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// Launch the PetscSolver
	solver->solve();

	solver->finalize();

	// Clean up
	MPI_Finalize();

	return;
}
