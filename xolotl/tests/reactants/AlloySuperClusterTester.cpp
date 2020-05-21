#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <AlloyCluster.h>
#include <AlloySuperCluster.h>
#include <AlloyClusterNetworkLoader.h>
#include <XolotlConfig.h>
#include <xolotlPerf/dummy/DummyHandlerRegistry.h>
#include <Constants.h>
#include <Options.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the AlloySuperCluster.
 */
BOOST_AUTO_TEST_SUITE(AlloySuperCluster_testSuite)

/**
 * This operation checks the ability of the AlloySuperCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
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

	// Create the loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(101);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& reactant = network->getAll(ReactantType::VoidSuper).begin()->second;

	// Check the type name
	BOOST_REQUIRE(ReactantType::VoidSuper == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for it
	int connectivityExpected[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the AlloySuperCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
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
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(101);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);
	// Add a grid point for the rates
	network->addGridPoints(1);
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Get the super cluster
	auto& cluster = network->getAll(ReactantType::VoidSuper).begin()->second;

	// Get one that it combines with (V1)
	auto secondCluster = (AlloyCluster *) network->get(Species::V, 1);
	// Set the temperature and concentration
	network->setTemperature(1000.0, 0);
	cluster->setConcentration(0.5);
	secondCluster->setConcentration(0.5);

	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(319615.775775, flux, 0.000001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the AlloySuperCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 638938, 246271, 111533, 42445.7, 0, 2.15429e+09,
			8.30345e+08, 3.76053e+08, 1.43113e+08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			-73.0464, 71.181, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 73.0464,
			-71.181, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
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
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(101);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Set the temperature in the network
	double temperature = 1000.0;
	network->setTemperature(temperature, 0);
	// Redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& cluster = network->getAll(ReactantType::VoidSuper).begin()->second;
	// Set the cluster concentration
	cluster->setConcentration(0.5);

	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 560U);

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
 * This operation checks the reaction radius for AlloySuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
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
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(101);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Check the reaction radius of the super cluster
	auto& cluster = network->getAll(ReactantType::VoidSuper).begin()->second;
	BOOST_REQUIRE_CLOSE(0.813023, cluster->getReactionRadius(), 0.001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
