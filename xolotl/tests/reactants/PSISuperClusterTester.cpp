#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include <PSISuperCluster.h>
#include <HDF5NetworkLoader.h>
#include <PSIHeCluster.h>
#include <PSIMixedCluster.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <Constants.h>
#include <Options.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the PSISuperCluster.
 */
BOOST_AUTO_TEST_SUITE(PSISuperCluster_testSuite)

/**
 * This operation checks the ability of the PSISuperCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(1);
	loader.setWidth(4, 0);
	loader.setWidth(1, 3);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Set the temperature in the network
	int networkSize = network->size();
	double temperature = 1000.0;
	network->setTemperature(temperature);
	network->computeRateConstants();
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Get the super cluster
	auto& reactant = network->getAll(ReactantType::PSISuper).begin()->second;

	// Check the type name
	BOOST_REQUIRE(ReactantType::PSISuper == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0,
			1, 1, 1, 1, 0 };
	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the ability of the PSISuperCluster to compute the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(1);
	loader.setWidth(4, 0);
	loader.setWidth(1, 3);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Set the temperature in the network
	int networkSize = network->size();
	double temperature = 1000.0;
	network->setTemperature(temperature);
	network->computeRateConstants();
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Get the super cluster
	auto& cluster = network->getAll(ReactantType::PSISuper).begin()->second;
	// Get one that it combines with (He)
	auto secondCluster = (PSICluster *) network->get(Species::He, 1);
	// Set the concentrations
	cluster->setConcentration(0.5);
	secondCluster->setConcentration(0.5);

	// Get and check the flux
	double flux = cluster->getTotalFlux();
	BOOST_REQUIRE_CLOSE(0.0, flux, 0.1);

	return;
}

/**
 * This operation checks the PSISuperCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(1);
	loader.setWidth(4, 0);
	loader.setWidth(1, 3);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Set the temperature in the network
	int networkSize = network->size();
	double temperature = 1000.0;
	network->setTemperature(temperature);
	network->computeRateConstants();
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Get the super cluster
	auto& cluster = network->getAll(ReactantType::PSISuper).begin()->second;

	// The vector of partial derivatives to compare with
	double knownPartials[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.24076e-19, 0.0, 0.0, 0.0 };
	// Set the concentration
	cluster->setConcentration(0.5);

	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives();

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 20U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	return;
}

/**
 * This operation checks the reaction radius for PSISuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(1);
	loader.setWidth(4, 0);
	loader.setWidth(1, 3);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Get the super cluster
	auto& cluster = network->getAll(ReactantType::PSISuper).begin()->second;

	// Check the radius
	BOOST_REQUIRE_CLOSE(0.137265, cluster->getReactionRadius(), 0.001);

	return;
}

/**
 * This operation checks the get concentration methods for PSISuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkGetConcentrations) {

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(1);
	loader.setWidth(4, 0);
	loader.setWidth(1, 3);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Get a super cluster
	auto& cluster = network->getAll(ReactantType::PSISuper).begin()->second;
	auto& superCluster =
			static_cast<PSISuperCluster&>(*cluster);

	// Set the concentration in the cluster
	superCluster.setZerothMoment(0.5);
	superCluster.setMoment(1.0, 0); // He
	superCluster.setMoment(1.0, 3); // V

	// Check the different concentrations
	BOOST_REQUIRE_CLOSE(0.5, superCluster.getTotalConcentration(), 0.001);
	BOOST_REQUIRE_CLOSE(4.5, superCluster.getTotalAtomConcentration(0), 0.001); // He
	BOOST_REQUIRE_CLOSE(0.5, superCluster.getTotalVacancyConcentration(), 0.001);
	BOOST_REQUIRE_CLOSE(0.5, superCluster.getIntegratedVConcentration(1), 0.001);

	return;
}

/**
 * This operation checks the boundary methods for PSISuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkBoundaries) {

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive_2D.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(1);
	loader.setWidth(4, 0);
	loader.setWidth(1, 3);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Get a super cluster
	auto& cluster = network->getAll(ReactantType::PSISuper).begin()->second;
	auto& superCluster =
			static_cast<PSISuperCluster&>(*cluster);

	// Check the different numbers
	BOOST_REQUIRE_CLOSE(1.0, superCluster.getNumV(), 0.001);
	BOOST_REQUIRE_CLOSE(1.0, superCluster.getNTot(), 0.001);
	auto bounds = superCluster.getBounds(0); // He
	BOOST_REQUIRE_EQUAL(9, *(bounds.begin()));
	BOOST_REQUIRE_EQUAL(10, *(bounds.end()));
	bounds = superCluster.getBounds(1); // D
	BOOST_REQUIRE_EQUAL(0, *(bounds.begin()));
	BOOST_REQUIRE_EQUAL(1, *(bounds.end()));
	bounds = superCluster.getBounds(2); // T
	BOOST_REQUIRE_EQUAL(0, *(bounds.begin()));
	BOOST_REQUIRE_EQUAL(1, *(bounds.end()));
	bounds = superCluster.getBounds(3); // V
	BOOST_REQUIRE_EQUAL(1, *(bounds.begin()));
	BOOST_REQUIRE_EQUAL(2, *(bounds.end()));

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
