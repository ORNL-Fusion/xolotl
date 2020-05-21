#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <FeCluster.h>
#include <FeSuperCluster.h>
#include <FeClusterNetworkLoader.h>
#include <FeHeCluster.h>
#include <FeHeVCluster.h>
#include <XolotlConfig.h>
#include <xolotlPerf/dummy/DummyHandlerRegistry.h>
#include <Constants.h>
#include <Options.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the FeSuperCluster.
 */
BOOST_AUTO_TEST_SUITE (FeSuperCluster_testSuite)

/**
 * This operation checks the ability of the FeSuperCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=6 0 0 6 1" << std::endl;
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
	// Initialize MPI for HDF5
	MPI_Init(&argc, &argv);

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the network loader
	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setVMin(4);
	loader.setHeWidth(2);
	loader.setVWidth(2);

	// Load the network
	auto network = loader.generate(opts);

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the first super cluster
	auto& reactant = network->getAll(ReactantType::FeSuper).begin()->second;

	// Check the type name
	BOOST_REQUIRE(ReactantType::FeSuper == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = { 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
			0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0,
			0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the ability of the FeSuperCluster to compute the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=6 0 0 6 1" << std::endl;
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

	// Create the network loader
	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setVMin(4);
	loader.setHeWidth(2);
	loader.setVWidth(2);
	// Load the network
	auto network = loader.generate(opts);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Set the temperature in the network
	double temperature = 1000.0;
	network->setTemperature(temperature, 0);
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& cluster = network->getAll(ReactantType::FeSuper).begin()->second;
	// Get one that it combines with (He)
	auto secondCluster = (FeCluster *) network->get(Species::He, 1);
	// Set the concentrations
	cluster->setConcentration(0.5);
	secondCluster->setConcentration(0.5);

	// Get and check the flux
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(43000201855.9, flux, 0.1);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the FeSuperCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=6 0 0 6 1" << std::endl;
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

	// Create the network loader
	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setVMin(4);
	loader.setHeWidth(2);
	loader.setVWidth(2);
	// Load the network
	auto network = loader.generate(opts);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Set the temperature in the network
	double temperature = 1000.0;
	network->setTemperature(temperature, 0);
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& cluster = network->getAll(ReactantType::FeSuper).begin()->second;

	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 2.37284e+09, 8.60004e+10, 0, 0, 0, 0, 0, 0, 0,
			5.15396e+07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			8.63906e-08, 0, 0, -0.00528717, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			8.63906e-08, 0, 0, 0, 0, 9.20495e-13, 0.00528717, 0 };
	// Set the concentration
	cluster->setConcentration(0.5);

	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 49U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the reaction radius for FeSuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=6 0 0 6 1" << std::endl;
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

	// Create the network loader
	FeClusterNetworkLoader loader = FeClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setVMin(4);
	loader.setHeWidth(2);
	loader.setVWidth(2);
	// Load the network
	auto network = loader.generate(opts);

	// Check the reaction connectivity of the super cluster
	auto& cluster = network->getAll(ReactantType::FeSuper).begin()->second;

	// Check the radius
	BOOST_REQUIRE_CLOSE(0.2492086, cluster->getReactionRadius(), 0.001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
