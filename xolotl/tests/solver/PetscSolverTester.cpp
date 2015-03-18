#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <string.h>
#include <PSIClusterNetworkLoader.h>
#include <PSIClusterReactionNetwork.h>
#include <PetscSolver.h>
#include <XolotlConfig.h>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>
#include <HDF5NetworkLoader.h>
#include <Options.h>
#include <PetscSolver1DHandler.h>
#include <PetscSolver2DHandler.h>
#include <PetscSolver3DHandler.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite configuration
 */BOOST_AUTO_TEST_SUITE (PetscSolverTester_testSuite)

/**
 * This operation checks the concentration of clusters after solving a test case
 * in 1D.
 */
BOOST_AUTO_TEST_CASE(checkPetscSolver1DHandler) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Local Declarations
	string sourceDir(XolotlSourceDirectory);

	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "vizHandler=dummy" << std::endl
			<< "petscArgs=-fieldsplit_0_pc_type redundant "
					"-ts_max_snes_failures 200 "
					"-pc_fieldsplit_detect_coupling "
					"-ts_adapt_dt_max 10 "
					"-pc_type fieldsplit "
					"-fieldsplit_1_pc_type sor "
					"-ts_final_time 1000 "
					"-ts_max_steps 5" << std::endl
			<< "startTemp=900" << std::endl
			<< "perfHandler=dummy" << std::endl
			<< "heFlux=4.0e5" << std::endl
			<< "material=W100" << std::endl
			<< "dimensions=1" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	argc = 1;
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the network loader
	std::shared_ptr<HDF5NetworkLoader> loader = std::make_shared<HDF5NetworkLoader>(
			make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Create the path to the network file
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string networkFilename = sourceDir + pathToFile;

	BOOST_TEST_MESSAGE(
			"PetscSolverTester Message: Network filename is: " << networkFilename);

	// Give the filename to the network loader
	loader->setFilename(networkFilename);

	// Create the solver
	std::shared_ptr<xolotlSolver::PetscSolver> solver = std::make_shared<
			xolotlSolver::PetscSolver>(make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Create the material factory
	auto materialFactory = xolotlFactory::IMaterialFactory::createMaterialFactory(opts.getMaterial(),
			opts.getDimensionNumber());

	// Initialize and get the temperature handler
	bool tempInitOK = xolotlFactory::initializeTempHandler(opts);
	auto tempHandler = xolotlFactory::getTemperatureHandler();

	// Set up our dummy performance and visualization infrastructures
    xolotlPerf::initialize(xolotlPerf::toPerfRegistryType("dummy"));
    xolotlFactory::initializeVizHandler(false);

    // Create a solver handler and initialize it
	auto solvHandler = std::make_shared<xolotlSolver::PetscSolver1DHandler> ();
	solvHandler->initializeHandlers(materialFactory, tempHandler, opts);

    // Set the solver command line to give the PETSc options and initialize it
    solver->setCommandLineOptions(opts.getPetscArgc(), opts.getPetscArgv());
	solver->initialize(solvHandler);

	// Give it the network loader
	solver->setNetworkLoader(loader);

	// Solve and finalize
	solver->solve();
	solver->finalize();

	// Check the concentrations left in the network
	auto network = solvHandler->getNetwork();
	double concs[network->getAll()->size()];
	network->fillConcentrationsArray(concs);

	// Check some concentrations
    BOOST_REQUIRE_CLOSE(concs[0], 3.9559e-11, 0.01);
    BOOST_REQUIRE_CLOSE(concs[1], 2.2999e-18, 0.01);
    BOOST_REQUIRE_CLOSE(concs[2], 4.1979e-26, 0.01);
    BOOST_REQUIRE_CLOSE(concs[7], -1.3576e-62, 0.01);
    BOOST_REQUIRE_CLOSE(concs[8], 0.0, 0.01);

    // Remove the created file
    std::string tempFile = "param.txt";
    std::remove(tempFile.c_str());
}

 /**
  * This operation checks the concentration of clusters after solving a test case
  * in 2D.
  */
 BOOST_AUTO_TEST_CASE(checkPetscSolver2DHandler) {
 	// Initialize MPI for HDF5
 	int argc = 0;
 	char **argv;
 	MPI_Init(&argc, &argv);

 	// Local Declarations
 	string sourceDir(XolotlSourceDirectory);

	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "vizHandler=dummy" << std::endl
			<< "petscArgs=-fieldsplit_0_pc_type redundant "
					"-ts_max_snes_failures 200 "
					"-pc_fieldsplit_detect_coupling "
					"-ts_adapt_dt_max 10 "
					"-pc_type fieldsplit "
					"-fieldsplit_1_pc_type sor "
					"-ts_final_time 1000 "
					"-ts_max_steps 5" << std::endl
			<< "startTemp=900" << std::endl
			<< "perfHandler=dummy" << std::endl
			<< "heFlux=4.0e5" << std::endl
			<< "material=W100" << std::endl
			<< "dimensions=2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	argc = 1;
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array

 	// Read the options
 	Options opts;
 	opts.readParams(argc, argv);

 	// Create the network loader
 	std::shared_ptr<HDF5NetworkLoader> loader = std::make_shared<HDF5NetworkLoader>(
 			make_shared<xolotlPerf::DummyHandlerRegistry>());

 	// Create the path to the network file
 	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
 	string networkFilename = sourceDir + pathToFile;

 	BOOST_TEST_MESSAGE(
 			"PetscSolverTester Message: Network filename is: " << networkFilename);

 	// Give the filename to the network loader
 	loader->setFilename(networkFilename);

 	// Create the solver
 	std::shared_ptr<xolotlSolver::PetscSolver> solver = std::make_shared<
 			xolotlSolver::PetscSolver>(make_shared<xolotlPerf::DummyHandlerRegistry>());

 	// Create the material factory
 	auto materialFactory = xolotlFactory::IMaterialFactory::createMaterialFactory(opts.getMaterial(),
 			opts.getDimensionNumber());

 	// Initialize and get the temperature handler
 	bool tempInitOK = xolotlFactory::initializeTempHandler(opts);
 	auto tempHandler = xolotlFactory::getTemperatureHandler();

 	// Set up our dummy performance and visualization infrastructures
     xolotlPerf::initialize(xolotlPerf::toPerfRegistryType("dummy"));
     xolotlFactory::initializeVizHandler(false);

     // Create a solver handler and initialize it
 	auto solvHandler = std::make_shared<xolotlSolver::PetscSolver2DHandler> ();
 	solvHandler->initializeHandlers(materialFactory, tempHandler, opts);

     // Set the solver command line to give the PETSc options and initialize it
     solver->setCommandLineOptions(opts.getPetscArgc(), opts.getPetscArgv());
 	solver->initialize(solvHandler);

 	// Give it the network loader
 	solver->setNetworkLoader(loader);

 	// Solve and finalize
 	solver->solve();
 	solver->finalize();

 	// Check the concentrations left in the network
 	auto network = solvHandler->getNetwork();
 	double concs[network->getAll()->size()];
 	network->fillConcentrationsArray(concs);

 	// Check some concentrations
     BOOST_REQUIRE_CLOSE(concs[0], 2.7841e-74, 0.01);
     BOOST_REQUIRE_CLOSE(concs[1], 1.6807e-146, 0.01);
     BOOST_REQUIRE_CLOSE(concs[6], 6.5151e-11, 0.01);
     BOOST_REQUIRE_CLOSE(concs[14], 0.0, 0.01);
     BOOST_REQUIRE_CLOSE(concs[23], 2.1984e-72, 0.01);

     // Remove the created file
     std::string tempFile = "param.txt";
     std::remove(tempFile.c_str());
 }

 /**
  * This operation checks the concentration of clusters after solving a test case
  * in 3D.
  */
 BOOST_AUTO_TEST_CASE(checkPetscSolver3DHandler) {
 	// Initialize MPI for HDF5
 	int argc = 0;
 	char **argv;
 	MPI_Init(&argc, &argv);

 	// Local Declarations
 	string sourceDir(XolotlSourceDirectory);

	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "vizHandler=dummy" << std::endl
			<< "petscArgs=-fieldsplit_0_pc_type redundant "
					"-ts_max_snes_failures 200 "
					"-pc_fieldsplit_detect_coupling "
					"-ts_adapt_dt_max 10 "
					"-pc_type fieldsplit "
					"-fieldsplit_1_pc_type sor "
					"-ts_final_time 1000 "
					"-ts_max_steps 5" << std::endl
			<< "startTemp=900" << std::endl
			<< "perfHandler=dummy" << std::endl
			<< "heFlux=4.0e5" << std::endl
			<< "material=W100" << std::endl
			<< "dimensions=3" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	argc = 1;
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array

 	// Read the options
 	Options opts;
 	opts.readParams(argc, argv);

 	// Create the network loader
 	std::shared_ptr<HDF5NetworkLoader> loader = std::make_shared<HDF5NetworkLoader>(
 			make_shared<xolotlPerf::DummyHandlerRegistry>());

 	// Create the path to the network file
 	string pathToFile("/tests/testfiles/tungsten_diminutive_3D.h5");
 	string networkFilename = sourceDir + pathToFile;

 	BOOST_TEST_MESSAGE(
 			"PetscSolverTester Message: Network filename is: " << networkFilename);

 	// Give the filename to the network loader
 	loader->setFilename(networkFilename);

 	// Create the solver
 	std::shared_ptr<xolotlSolver::PetscSolver> solver = std::make_shared<
 			xolotlSolver::PetscSolver>(make_shared<xolotlPerf::DummyHandlerRegistry>());

 	// Create the material factory
 	auto materialFactory = xolotlFactory::IMaterialFactory::createMaterialFactory(opts.getMaterial(),
 			opts.getDimensionNumber());

 	// Initialize and get the temperature handler
 	bool tempInitOK = xolotlFactory::initializeTempHandler(opts);
 	auto tempHandler = xolotlFactory::getTemperatureHandler();

 	// Set up our dummy performance and visualization infrastructures
     xolotlPerf::initialize(xolotlPerf::toPerfRegistryType("dummy"));
     xolotlFactory::initializeVizHandler(false);

     // Create a solver handler and initialize it
 	auto solvHandler = std::make_shared<xolotlSolver::PetscSolver3DHandler> ();
 	solvHandler->initializeHandlers(materialFactory, tempHandler, opts);

     // Set the solver command line to give the PETSc options and initialize it
     solver->setCommandLineOptions(opts.getPetscArgc(), opts.getPetscArgv());
 	solver->initialize(solvHandler);

 	// Give it the network loader
 	solver->setNetworkLoader(loader);

 	// Solve and finalize
 	solver->solve();
 	solver->finalize();

 	// Check the concentrations left in the network
 	auto network = solvHandler->getNetwork();
 	double concs[network->getAll()->size()];
 	network->fillConcentrationsArray(concs);

 	// Check some concentrations
     BOOST_REQUIRE_CLOSE(concs[0], -7.4983e-79, 0.01);
     BOOST_REQUIRE_CLOSE(concs[6], 8.8196e-11, 0.01);
     BOOST_REQUIRE_CLOSE(concs[14], 0.0, 0.01);
     BOOST_REQUIRE_CLOSE(concs[15], -3.0887e-157, 0.01);
     BOOST_REQUIRE_CLOSE(concs[16], -5.5119e-144, 0.01);

     // Remove the created file
     std::string tempFile = "param.txt";
     std::remove(tempFile.c_str());
 }

BOOST_AUTO_TEST_SUITE_END()
