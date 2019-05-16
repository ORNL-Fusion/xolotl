#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <NECluster.h>
#include <NESuperCluster.h>
#include <NEClusterNetworkLoader.h>
#include <NEXeCluster.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <Constants.h>
#include <Options.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the NESuperCluster.
 */
BOOST_AUTO_TEST_SUITE(NESuperCluster_testSuite)

/**
 * This operation checks the ability of the NESuperCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
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
	// Initialize MPI for HDF5
	MPI_Init(&argc, &argv);

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

	// Redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& reactant = network->getAll(ReactantType::NESuper).begin()->second;

	// Check the type name
	BOOST_REQUIRE(ReactantType::NESuper == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Xe
	int connectivityExpected[] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the NESuperCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
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
	// Add a grid point for the rates
	network->addGridPoints(1);
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Get the super cluster
	auto& cluster = network->getAll(ReactantType::NESuper).begin()->second;

	// Get one that it combines with (Xe1)
	auto secondCluster = (NECluster *) network->get(Species::Xe, 1);
	// Set the temperature and concentration
	network->setTemperature(1000.0, 0);
	cluster->setConcentration(0.5);
	secondCluster->setConcentration(0.5);

	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(0.0, flux, 0.000001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the NESuperCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 0.0540719, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			-1.18392e-36, 1.20469e-36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			1.18392e-36, 0, 0 };

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
	// Set a fission rate for the diffusion to work
	network->setFissionRate(8.0e-9);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Set the temperature in the network
	double temperature = 1000.0;
	network->setTemperature(temperature, 0);
	// Redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& cluster = network->getAll(ReactantType::NESuper).begin()->second;
	// Set the cluster concentration
	cluster->setConcentration(0.5);

	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 70U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.001);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the reaction radius for NESuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
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

	// Check the reaction radius of the super cluster
	auto& cluster = network->getAll(ReactantType::NESuper).begin()->second;
	BOOST_REQUIRE_CLOSE(1.3135906803, cluster->getReactionRadius(), 0.001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
