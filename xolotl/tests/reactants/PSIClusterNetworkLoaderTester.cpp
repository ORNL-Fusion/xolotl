#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <DummyHandlerRegistry.h>
#include <PSIClusterReactionNetwork.h>
#include <Options.h>
#include "tests/utils/MPIFixture.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

// Initialize MPI before running any tests; finalize it running all tests.
BOOST_GLOBAL_FIXTURE(MPIFixture);

/**
 * This suite is responsible for testing the PSIClusterNetworkLoader. It
 * creates a string stream that contains each of the available PSICluster types
 * and checks that the loader returns a list with each type in it.
 */
BOOST_AUTO_TEST_SUITE(PSIClusterNetworkLoader_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkLoading) {
	// Local Declarations
	shared_ptr<stringstream> networkStream(
			new stringstream(stringstream::in | stringstream::out));
	string singleHeString = "1 0 0 0 0 6.15 0.999 1.34\n";
	string singleVString = "0 0 0 50 0 3.6 0.888 2.345\n";
	string singleIString = "0 0 0 0 1 5.0 0.7777 3.456\n";
	string mixedString = "1 0 0 50 0 2.49 6.789 4.5678\n";
	// This string is bad because it is one value short
	string badString = "1 2 3\n";
	bool caughtFlag = false;
	PSIClusterNetworkLoader loader = PSIClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Load the network stream. This simulates a file with single He, single
	// V, single I and one mixed-species cluster. They are mixed up here to test
	// the ability of the loader to order them.
	*networkStream << singleVString << mixedString << singleHeString
			<< singleIString;

	// Setup the Loader
	loader.setInputstream(networkStream);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);
	auto psiNetwork = (PSIClusterReactionNetwork*) network.get();

	// Check the properties
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::He) == 1);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::V) == 50);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::I) == 1);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::PSIMixed) == 51);

	// Check the reactants - He first
	auto heCluster = (PSICluster *) network->get(Species::He, 1);
	BOOST_REQUIRE(heCluster->getSize() == 1);
	double formationEnergy = heCluster->getFormationEnergy();
	BOOST_REQUIRE_CLOSE(formationEnergy, 6.15, 0.001);
	// V
	auto vCluster = (PSICluster *) network->get(Species::V, 50);
	BOOST_REQUIRE(vCluster->getSize() == 50);
	formationEnergy = vCluster->getFormationEnergy();
	BOOST_REQUIRE_CLOSE(formationEnergy, 3.6, 0.001);
	BOOST_REQUIRE_CLOSE(vCluster->getMigrationEnergy(), 0.888, 0.001);
	BOOST_REQUIRE_CLOSE(vCluster->getDiffusionFactor(), 2.345, 0.001);
	// I
	auto iCluster = (PSICluster *) network->get(Species::I, 1);
	BOOST_REQUIRE(iCluster->getSize() == 1);
	formationEnergy = iCluster->getFormationEnergy();
	BOOST_REQUIRE_CLOSE(formationEnergy, 5.0, 0.001);
	BOOST_REQUIRE_CLOSE(iCluster->getMigrationEnergy(), 0.7777, 0.0001);
	BOOST_REQUIRE_CLOSE(iCluster->getDiffusionFactor(), 3.456, 0.001);
	// HeV
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 1;
	composition[toCompIdx(Species::V)] = 50;
	composition[toCompIdx(Species::I)] = 0;
	auto heVCluster = (PSICluster *) network->get(ReactantType::PSIMixed,
			composition);
	BOOST_REQUIRE(heVCluster->getSize() == 51);
	formationEnergy = heVCluster->getFormationEnergy();
	BOOST_REQUIRE_CLOSE(formationEnergy, 2.49, 0.001);
	BOOST_REQUIRE_CLOSE(heVCluster->getMigrationEnergy(), 6.789, 0.001);
	BOOST_REQUIRE_CLOSE(heVCluster->getDiffusionFactor(), 4.5678, 0.0001);

	// Reload the network stream with the bad string
	(*networkStream).clear();
	*networkStream << badString;
	// Make sure the exception is caught when loading the bad string
	try {
		loader.load(opts);
	} catch (const string& /* error */) {
		// Do nothing but flip the flag
		caughtFlag = true;
	}
	BOOST_REQUIRE(caughtFlag);
}

/**
 * Method checking the generation of the network.
 */
BOOST_AUTO_TEST_CASE(checkGenerate) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 5 3" << std::endl << "grid=100 0.5"
			<< std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the loader
	PSIClusterNetworkLoader loader = PSIClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 104);

	// Check the properties
	auto psiNetwork = (PSIClusterReactionNetwork*) network.get();
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::He) == 8);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::V) == 5);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::I) == 3);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::PSIMixed) == 32);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

BOOST_AUTO_TEST_SUITE_END()
