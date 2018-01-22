#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <NECluster.h>
#include <NESuperCluster.h>
#include <NEClusterNetworkLoader.h>
#include <XeCluster.h>
#include <XolotlConfig.h>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>
#include <Constants.h>

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
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/fuel_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setXeMin(2);
	loader.setWidth(2);

	// Load the network
	auto network = loader.load();

	// Set the temperature in the network
	int networkSize = network->size();
	auto allReactants = network->getAll();
	double temperature = 1000.0;
	for (int i = 0; i < networkSize; i++) {
		// This part will set the temperature in each reactant
		// and recompute the diffusion coefficient
		allReactants->at(i)->setTemperature(temperature);
	}
	network->computeRateConstants();
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto reactant = network->getAll(NESuperType).at(0);

	// Check the type name
	BOOST_REQUIRE_EQUAL(NESuperType, reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Xe
	int connectivityExpected[] = { 1, 1, 1, 0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the NESuperCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {

	// Create the network loader
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/fuel_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setXeMin(2);
	loader.setWidth(2);

	// Load the network
	auto network = loader.load();

	// Check the reaction connectivity of the super cluster
	auto cluster = (NECluster *) network->getAll(NESuperType).at(0);

	// Get one that it combines with (Xe1)
	auto secondCluster = (NECluster *) network->get(xeType, 1);
	// Set the temperature and concentration
	network->setTemperature(1000.0);
	cluster->setConcentration(0.5);
	secondCluster->setConcentration(0.5);

	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux();
	BOOST_TEST_MESSAGE(
			"XeClusterTester Message: \n" << "Total Flux is " << flux << "\n" << "   -Production Flux: " << cluster->getProductionFlux() << "\n" << "   -Combination Flux: " << cluster->getCombinationFlux() << "\n" << "   -Dissociation Flux: " << cluster->getDissociationFlux() << "\n" << "   -Emission Flux: " << cluster->getEmissionFlux() << "\n");

	BOOST_REQUIRE_CLOSE(0.00942477796, flux, 0.000001);

	return;
}

/**
 * This operation checks the NESuperCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 0.0, -752.45682, 752.45682, 0.0 };

	// Create the network loader
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/fuel_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setXeMin(2);
	loader.setWidth(2);

	// Load the network
	auto network = loader.load();

	// Check the reaction connectivity of the super cluster
	auto cluster = (NECluster *) network->getAll(NESuperType).at(0);

	// Set the temperature in the network
	int networkSize = network->size();
	auto allReactants = network->getAll();
	double temperature = 1000.0;
	network->setTemperature(temperature);
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();
	// Set the cluster concentration
	cluster->setConcentration(0.5);
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives();

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 4U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.001);
	}

	return;
}

/**
 * This operation checks the reaction radius for NESuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	// Create the network loader
	NEClusterNetworkLoader loader = NEClusterNetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/fuel_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setXeMin(2);
	loader.setWidth(2);

	// Load the network
	auto network = loader.load();

	// Check the reaction radius of the super cluster
	auto cluster = network->getAll(NESuperType).at(0);
	BOOST_REQUIRE_CLOSE(0.3869446, cluster->getReactionRadius(), 0.001);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
