#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <NEClusterReactionNetwork.h>
#include <NEClusterNetworkLoader.h>
#include <NECluster.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the NEClusterNetworkLoader.
 */
BOOST_AUTO_TEST_SUITE(NEClusterNetworkLoader_testSuite)

/**
 * Method checking the loading of the network from the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkLoad) {
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

	// Load the network
	auto network = loader.load();
	auto neNetwork = std::dynamic_pointer_cast<NEClusterReactionNetwork>(
			network);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 3);

	// Check the properties
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxXeClusterSize(), 3);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxVClusterSize(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxIClusterSize(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxXeVClusterSize(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumXeClusters(), 3);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumVClusters(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumIClusters(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumSuperClusters(), 0);

	// Get all the reactants
	auto reactants = network->getAll();

	// Get the first one of the network
	auto reactant = (NECluster *) reactants->at(0);
	// Check the composition
	auto composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["Xe"], 1);
	BOOST_REQUIRE_EQUAL(composition["V"], 0);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
	// Check the formation energy
	auto formationEnergy = reactant->getFormationEnergy();
	BOOST_REQUIRE_EQUAL(formationEnergy, 7.0);
	// Check the migration energy
	auto migrationEnergy = reactant->getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.0);
	// Check the diffusion factor
	auto diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 5.0e-3);

	// Get the last reactant of the network
	reactant = (NECluster *) reactants->at(2);
	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["Xe"], 3);
	BOOST_REQUIRE_EQUAL(composition["V"], 0);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
	// Check the formation energy
	formationEnergy = reactant->getFormationEnergy();
	BOOST_REQUIRE_EQUAL(formationEnergy, 17.15);
	// Check the diffusion factor
	diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_CLOSE(diffusionFactor, 0.0, 1.0e-16);

	return;
}

/**
 * Method checking the loading of the network from the HDF5 file and
 * the apply sectional method.
 */
BOOST_AUTO_TEST_CASE(checkApplySectional) {

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

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 2);

	// Get the dof of the network
	int dof = network->getDOF();
	// Check the value
	BOOST_REQUIRE_EQUAL(dof, 3);

	// Check the properties
	auto neNetwork = std::dynamic_pointer_cast<NEClusterReactionNetwork>(
			network);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxXeClusterSize(), 3);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxVClusterSize(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxIClusterSize(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxXeVClusterSize(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumXeClusters(), 1);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumVClusters(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumIClusters(), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getNumSuperClusters(), 1);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
