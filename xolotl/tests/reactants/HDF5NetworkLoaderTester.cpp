#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <PSIClusterReactionNetwork.h>
#include <DummyHandlerRegistry.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>
#include <Options.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HDF5NetworkLoader.
 */
BOOST_AUTO_TEST_SUITE(HDF5NetworkLoader_testSuite)

/**
 * Method checking the loading of the network from the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkLoad) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 9);

	// Check the properties
	auto psiNetwork = (PSIClusterReactionNetwork*) network.get();
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::He) == 8);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::V) == 1);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::I) == 0);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::PSIMixed) == 0);

	// Get all the reactants
	auto& reactants = network->getAll();

	// Get the first one of the network
	IReactant& reactant = reactants.at(0);
	// Check the composition
	auto composition = reactant.getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 1);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);
	// Check the formation energy
	auto formationEnergy = reactant.getFormationEnergy();
	BOOST_REQUIRE_EQUAL(formationEnergy, 6.15);
	// Check the migration energy
	auto migrationEnergy = reactant.getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.13);
	// Check the diffusion factor
	auto diffusionFactor = reactant.getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 2.9e+10);

	// Get the last reactant of the network
	IReactant& reactantBis = reactants.at(8);
	// Check the composition
	composition = reactantBis.getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 1);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);
	// Check the formation energy
	formationEnergy = reactantBis.getFormationEnergy();
	BOOST_REQUIRE_EQUAL(formationEnergy, 3.6);
	// Check the migration energy
	migrationEnergy = reactantBis.getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 1.3);
	// Check the diffusion factor
	diffusionFactor = reactantBis.getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 1.8e+12);

	return;
}

/**
 * Method checking the loading of the network from the HDF5 file and
 * the apply sectional method.
 */
BOOST_AUTO_TEST_CASE(checkApplySectional) {

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);
	// Set grouping parameters
	loader.setVMin(28);
	loader.setWidth(4, 0);
	loader.setWidth(2, 3);

	// Create default options
	Options opts;

	// Load the network
	auto network = loader.load(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 1870);

	// Get the dof of the network
	int dof = network->getDOF();
	// Check the value
	BOOST_REQUIRE_EQUAL(dof, 1933);

	// Check the properties
	auto psiNetwork = (PSIClusterReactionNetwork*) network.get();

	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::He) == 8);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::V) == 29);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::I) == 6);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::PSIMixed) == 137);
	BOOST_REQUIRE(psiNetwork->getMaxClusterSize(ReactantType::PSISuper) == 145);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
