#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <Options.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite configuration
 */
BOOST_AUTO_TEST_SUITE (TungstenIntegrationTester_testSuite)

/**
 * This operation checks the fluxs from the reactant as best as is possible
 * given that it requires external data.
 */
BOOST_AUTO_TEST_CASE(checkGetReactantFluxesAndParials) {
	// Local Declarations
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten.txt");
	string networkFilename = sourceDir + pathToFile;

	BOOST_TEST_MESSAGE(
			"TungstenIntegrationTester Message: Network filename is: " << networkFilename);

	// Load the input file from the master task
	shared_ptr<istream> networkStream = make_shared<ifstream>(networkFilename);

	// Create a network loader and set the istream on every MPI task
	shared_ptr<PSIClusterNetworkLoader> networkLoader = make_shared<
			PSIClusterNetworkLoader>(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	networkLoader->setInputstream(networkStream);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = networkLoader->load(opts);
	// Add rates
	network->addGridPoints(1);

	BOOST_TEST_MESSAGE("TungstenIntegrationTester Message: Network loaded");

	// Get all the reactants
	auto& allReactants = network->getAll();
	// Get the network size
	const int size = network->size();
	// Set the temperature
	double temperature = 1000.0;
	network->setTemperature(temperature, 0);

	// Initialize all the concentrations to 0.001;
	for (int i = 0; i < size; ++i) {
		IReactant& reactant = allReactants.at(i);
		reactant.setConcentration(0.001);
	}

	BOOST_TEST_MESSAGE(
			"TungstenIntegrationTester Message: " << "Size of the network is: " << size);

	// Create an array to test the second partial derivative routine.
	auto secondPartials = std::vector<double>(size, 0.0);

	BOOST_TEST_MESSAGE("Check partial derivatives.");
	for (int i = 0; i < size; ++i) {
		IReactant& reactant = allReactants.at(i);
		// Get the partials using method 1
		auto partials = reactant.getPartialDerivatives(0);
		// Get the partials using method 2
		reactant.getPartialDerivatives(secondPartials, 0);
		// Compare the two arrays of partial derivatives
		for (int j = 0; j < size; ++j) {
			BOOST_REQUIRE_CLOSE(partials[j], secondPartials[j], 1.0);
		}
		// Zero the partials array
		std::fill(secondPartials.begin(), secondPartials.end(), 0.0);
	}

	return;
}

/**
 * This operation checks the partial derivatives for the production
 * He_1 + He_1 --> He_2
 * and the dissociation
 * He_2 --> He_1 + He_1
 */
BOOST_AUTO_TEST_CASE(checkSingleReaction) {
	// Local Declarations
	string sourceDir(XolotlSourceDirectory);
	// This file contains specific values to obtain round numbers for the partial derivatives
	string pathToFile("/tests/testfiles/single_reaction.txt");
	string networkFilename = sourceDir + pathToFile;

	BOOST_TEST_MESSAGE(
			"TungstenIntegrationTester Message: Network filename is: " << networkFilename);

	// Load the input file from the master task
	shared_ptr<istream> networkStream = make_shared<ifstream>(networkFilename);

	// Create a network loader and set the istream on every MPI task
	shared_ptr<PSIClusterNetworkLoader> networkLoader = make_shared<
			PSIClusterNetworkLoader>(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	networkLoader->setInputstream(networkStream);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = networkLoader->load(opts);
	// Add rates
	network->addGridPoints(1);

	BOOST_TEST_MESSAGE("TungstenIntegrationTester Message: Network loaded");

	// Get all the reactants
	auto& allReactants = network->getAll();
	// Get the network size
	const int size = network->size();
	// Set the temperature
	double temperature = 1000.0;
	// Initialize the rate constants
	network->setTemperature(temperature, 0);

	// Initialize all the concentrations to 0.001;
	for (int i = 0; i < size; ++i) {
		IReactant& reactant = allReactants.at(i);
		reactant.setConcentration(0.001);
	}

	// Get He_1
	IReactant& reactant = allReactants.at(0);
	// Its partial derivatives
	auto partials = reactant.getPartialDerivatives(0);

	// Check the values of the partial derivatives
	BOOST_REQUIRE_CLOSE(partials[0], -2.0, 0.1);
	BOOST_REQUIRE_CLOSE(partials[1], 1.0, 0.1);

	// Get He_2
	IReactant& reactantBis = allReactants.at(1);
	// Its partial derivatives
	partials = reactantBis.getPartialDerivatives(0);

	// Check the values of the partial derivatives
	BOOST_REQUIRE_CLOSE(partials[0], 1.0, 0.1);
	BOOST_REQUIRE_CLOSE(partials[1], -0.5, 0.1);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
