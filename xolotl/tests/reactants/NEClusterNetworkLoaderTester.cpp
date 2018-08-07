#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <NEClusterReactionNetwork.h>
#include <NEClusterNetworkLoader.h>
#include <NECluster.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>
#include <Options.h>
#include "tests/utils/MPIFixture.h"

using namespace std;
using namespace xolotlCore;

// Initialize MPI before running any tests; finalize it running all tests.
BOOST_GLOBAL_FIXTURE(MPIFixture);

/**
 * This suite is responsible for testing the NEClusterNetworkLoader.
 */
BOOST_AUTO_TEST_SUITE(NEClusterNetworkLoader_testSuite)

/**
 * Method checking the loading of the network from the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkLoad) {

	// Create the network loader
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/fuel_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);
	auto neNetwork = (NEClusterReactionNetwork*) network.get();

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 0);

	// Check the properties
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::Xe), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::XeV), 0);

//	// Get all the reactants
//	auto reactants = network->getAll();
//
//	// Get the first one of the network
//	IReactant& reactant = reactants.at(0);
//	// Check the composition
//	auto composition = reactant.getComposition();
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Xe)], 1);
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
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::Xe)], 3);
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
//	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);
//	// Check the formation energy
//	formationEnergy = reactant2.getFormationEnergy();
//	BOOST_REQUIRE_EQUAL(formationEnergy, 17.15);
//	// Check the diffusion factor
//	diffusionFactor = reactant2.getDiffusionFactor();
//	BOOST_REQUIRE_CLOSE(diffusionFactor, 0.0, 1.0e-16);

	return;
}

/**
 * Method checking the generation of the network.
 */
BOOST_AUTO_TEST_CASE(checkGenerate) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=100" << std::endl << "grid=100 0.5" << std::endl;
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
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 100);

	// Check the properties
	auto neNetwork = (NEClusterReactionNetwork*) network.get();
	// Check the properties
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::Xe), 100);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::XeV), 0);

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
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=100" << std::endl << "grid=100 0.5" << std::endl;
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
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setXeMin(2);
	loader.setWidth(2);

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 51);

	// Get the dof of the network
	int dof = network->getDOF();
	// Check the value
	BOOST_REQUIRE_EQUAL(dof, 102);

	// Check the properties
	auto neNetwork = (NEClusterReactionNetwork*) network.get();
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::Xe), 100);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::XeV), 0);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

BOOST_AUTO_TEST_SUITE_END()
