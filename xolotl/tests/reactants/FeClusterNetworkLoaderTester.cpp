#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <FeClusterReactionNetwork.h>
#include <FeClusterNetworkLoader.h>
#include <FeCluster.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>
#include <Options.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the FeClusterNetworkLoader.
 */
BOOST_AUTO_TEST_SUITE(FeClusterNetworkLoader_testSuite)

///**
// * Method checking the loading of the network from the HDF5 file.
// */
//BOOST_AUTO_TEST_CASE(checkLoad) {
//	// Initialize MPI for HDF5
//	int argc = 0;
//	char **argv;
//	MPI_Init(&argc, &argv);
//
//	// Create the network loader
//	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
//			make_shared<xolotlPerf::DummyHandlerRegistry>());
//	// Define the filename to load the network from
//	string sourceDir(XolotlSourceDirectory);
//	string pathToFile("/tests/testfiles/fuel_diminutive.h5");
//	string filename = sourceDir + pathToFile;
//	// Give the filename to the network loader
//	loader.setFilename(filename);
//
//	// Create the options needed to load the network
//	Options opts;
//	// Load the network
//	auto network = loader.load(opts);
//	auto neNetwork = (FeClusterReactionNetwork*) network.get();
//
//	// Get the size of the network
//	int networkSize = network->size();
//	// Check the value
//	BOOST_REQUIRE_EQUAL(networkSize, 3);
//
//	// Check the properties
//	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::He), 3);
//	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 0);
//	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 0);
//	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::HeV), 0);
//
//	// Get all the reactants
//	auto reactants = network->getAll();
//
//	// Get the first one of the network
//	IReactant& reactant = reactants.at(0);
//	// Check the composition
//	auto composition = reactant.getComposition();
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 1);
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);
//	// Check the formation energy
//	auto formationEnergy = reactant.getFormationEnergy();
//	BOOST_REQUIRE_EQUAL(formationEnergy, 7.0);
//	// Check the migration energy
//	auto migrationEnergy = reactant.getMigrationEnergy();
//	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.0);
//	// Check the diffusion factor
//	auto diffusionFactor = reactant.getDiffusionFactor();
//	BOOST_REQUIRE_EQUAL(diffusionFactor, 5.0e-3);
//
//	// Get the last reactant of the network
//	IReactant& reactant2 = reactants.at(2);
//	// Check the composition
//	composition = reactant2.getComposition();
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 3);
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);
//	// Check the formation energy
//	formationEnergy = reactant2.getFormationEnergy();
//	BOOST_REQUIRE_EQUAL(formationEnergy, 17.15);
//	// Check the diffusion factor
//	diffusionFactor = reactant2.getDiffusionFactor();
//	BOOST_REQUIRE_CLOSE(diffusionFactor, 0.0, 1.0e-16);
//
//	return;
//}

/**
 * Method checking the generation of the network.
 */
BOOST_AUTO_TEST_CASE(checkGenerate) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=5 5 1" << std::endl << "grid=100 0.5" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argv);

	// Create the loader
	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 39);

	// Check the properties
	auto neNetwork = (FeClusterReactionNetwork*) network.get();
	// Check the properties
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::He), 8);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 5);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 1);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::HeV), 10);

	// Get all the reactants
	auto reactants = network->getAll();

	// Get the first one of the network
	IReactant& reactant = reactants.at(0);
	// Check the composition
	auto composition = reactant.getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 1);
	// Check the formation energy
	auto formationEnergy = reactant.getFormationEnergy();
	BOOST_REQUIRE_EQUAL(formationEnergy, 0.0);
	// Check the migration energy
	auto migrationEnergy = reactant.getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.34);
	// Check the diffusion factor
	auto diffusionFactor = reactant.getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 100000000000);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * Method checking the loading of the network from the HDF5 file and
 * the apply sectional method.
 */
BOOST_AUTO_TEST_CASE(checkApplySectional) {

	// Create the network loader
	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setVMin(4);
	loader.setHeWidth(2);
	loader.setVWidth(2);

	// Create the options needed to load the network
	Options opts;
	opts.setMaxV(6);
	opts.setMaxImpurity(6);
	opts.setMaxI(1);
	// Load the network
	auto network = loader.generate(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 32);

	// Get the dof of the network
	int dof = network->getDOF();
	// Check the value
	BOOST_REQUIRE_EQUAL(dof, 49);

	// Check the properties
	auto neNetwork = (FeClusterReactionNetwork*) network.get();
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::He), 8);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 6);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 1);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::HeV), 12);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
