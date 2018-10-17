#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <AlloyClusterReactionNetwork.h>
#include <AlloyClusterNetworkLoader.h>
#include <AlloyCluster.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>
#include <Options.h>
#include "tests/utils/MPIFixture.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

// Initialize MPI before running any tests; finalize it running all tests.
BOOST_GLOBAL_FIXTURE(MPIFixture);

/**
 * This suite is responsible for testing the AlloyClusterNetworkLoader.
 */
BOOST_AUTO_TEST_SUITE(AlloyClusterNetworkLoader_testSuite)

/**
 * Method checking the loading of the network from the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkLoad) {

	// Create the network loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/alloy_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);
	auto alloyNetwork = (AlloyClusterReactionNetwork*) network.get();

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 1115);

	// Check the properties
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::V), 5);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::I), 4);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Void), 1000);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Frank), 1000);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Faulted), 1000);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Perfect), 1000);

	// Get all the reactants
	auto reactants = network->getAll();

	// Get the first one of the network
	IReactant& reactant = reactants.at(0);
	// Check the composition
	auto composition = reactant.getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 1);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Void)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Frank)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Faulted)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Perfect)], 0);
	// Check the formation energy
	auto formationEnergy = reactant.getFormationEnergy();
	BOOST_REQUIRE_EQUAL(formationEnergy, 1.5);
	// Check the migration energy
	auto migrationEnergy = reactant.getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 1.2);
	// Check the diffusion factor
	auto diffusionFactor = reactant.getDiffusionFactor();
	BOOST_REQUIRE_CLOSE(diffusionFactor, 1.08e11, 1.0e-5);

	// Get the last reactant of the network
	IReactant& reactant2 = reactants.at(6);
	// Check the composition
	composition = reactant2.getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 2);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Void)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Frank)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Faulted)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Perfect)], 0);
	// Check the formation energy
	formationEnergy = reactant2.getFormationEnergy();
	BOOST_REQUIRE_CLOSE(formationEnergy, 6.0559, 0.001);
	// Check the diffusion factor
	diffusionFactor = reactant2.getDiffusionFactor();
	BOOST_REQUIRE_CLOSE(diffusionFactor, 5.4e10, 1.0e-5);
	// Check the migration energy
	migrationEnergy = reactant2.getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.5);

	return;
}

/**
 * Method checking the generation of the network.
 */
BOOST_AUTO_TEST_CASE(checkGenerate) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=100 0 0 5 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 0;
	char **argv;
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argv);

	// Create the loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 391);

	// Check the properties
	auto alloyNetwork = (AlloyClusterReactionNetwork*) network.get();
	// Check the properties
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::V), 5);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::I), 4);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Void), 100);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Frank), 100);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Faulted), 100);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Perfect), 100);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

BOOST_AUTO_TEST_SUITE_END()
