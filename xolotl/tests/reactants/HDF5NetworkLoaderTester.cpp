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
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = loader.load();

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 5);

	// Get the properties
	auto props = network->getProperties();

	// Check the properties
	BOOST_REQUIRE_EQUAL(strtol(props["maxHeClusterSize"].c_str(),NULL,10), 2);
	BOOST_REQUIRE_EQUAL(strtol(props["maxVClusterSize"].c_str(),NULL,10), 1);
	BOOST_REQUIRE_EQUAL(strtol(props["maxIClusterSize"].c_str(),NULL,10), 1);
	BOOST_REQUIRE_EQUAL(strtol(props["maxHeVClusterSize"].c_str(),NULL,10), 3);
	BOOST_REQUIRE_EQUAL(strtol(props["numHeClusters"].c_str(),NULL,10), 2);
	BOOST_REQUIRE_EQUAL(strtol(props["numVClusters"].c_str(),NULL,10), 1);
	BOOST_REQUIRE_EQUAL(strtol(props["numIClusters"].c_str(),NULL,10), 1);

	// Get all the reactants
	auto reactants = network->getAll();

	// Get the first one of the network
	auto reactant = (PSICluster *) reactants->at(0);
	// Check the composition
	auto composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["He"], 0);
	BOOST_REQUIRE_EQUAL(composition["V"], 1);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
	// Check the binding energies
	auto bindingEnergies = reactant->getBindingEnergies();
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(0), numeric_limits<double>::infinity()); // Helium binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(1), numeric_limits<double>::infinity()); // Vacancy binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(2), numeric_limits<double>::infinity()); // Interstitial binding energy
	// Check the migration energy
	auto migrationEnergy = reactant->getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 1.66);
	// Check the diffusion factor
	auto diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 2.41e+11);

	// Get the second reactant of the network
	reactant = (PSICluster *) reactants->at(1);
	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["He"], 2);
	BOOST_REQUIRE_EQUAL(composition["V"], 0);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
	// Check the binding energies
	bindingEnergies = reactant->getBindingEnergies();
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(0), 0.864); // Helium binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(1), numeric_limits<double>::infinity()); // Vacancy binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(2), numeric_limits<double>::infinity()); // Interstitial binding energy
	// Check the migration energy
	migrationEnergy = reactant->getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.2);
	// Check the diffusion factor
	diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 3.24e+10);

	// Get the third reactant of the network
	reactant = (PSICluster *) reactants->at(2);// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["He"], 1);
	BOOST_REQUIRE_EQUAL(composition["V"], 0);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
	// Check the binding energies
	bindingEnergies = reactant->getBindingEnergies();
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(0), numeric_limits<double>::infinity()); // Helium binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(1), numeric_limits<double>::infinity()); // Vacancy binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(2), numeric_limits<double>::infinity()); // Interstitial binding energy
	// Check the migration energy
	migrationEnergy = reactant->getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.13);
	// Check the diffusion factor
	diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 2.95e+10);

	// Get the fourth reactant of the network
	reactant = (PSICluster *) reactants->at(3);
	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["He"], 0);
	BOOST_REQUIRE_EQUAL(composition["V"], 0);
	BOOST_REQUIRE_EQUAL(composition["I"], 1);
	// Check the binding energies
	bindingEnergies = reactant->getBindingEnergies();
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(0), numeric_limits<double>::infinity()); // Helium binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(1), numeric_limits<double>::infinity()); // Vacancy binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(2), numeric_limits<double>::infinity()); // Interstitial binding energy
	// Check the migration energy
	migrationEnergy = reactant->getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, 0.013);
	// Check the diffusion factor
	diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 2.13e+10);

	// Get the last and fifth reactant of the network
	reactant = (PSICluster *) reactants->at(4);
	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition["He"], 2);
	BOOST_REQUIRE_EQUAL(composition["V"], 1);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
	// Check the binding energies
	bindingEnergies = reactant->getBindingEnergies();
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(0), 3.02); // Helium binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(1), 7.25); // Vacancy binding energy
	BOOST_REQUIRE_EQUAL(bindingEnergies.at(2), numeric_limits<double>::infinity()); // Interstitial binding energy
	// Check the migration energy
	migrationEnergy = reactant->getMigrationEnergy();
	BOOST_REQUIRE_EQUAL(migrationEnergy, numeric_limits<double>::infinity());
	// Check the diffusion factor
	diffusionFactor = reactant->getDiffusionFactor();
	BOOST_REQUIRE_EQUAL(diffusionFactor, 0.0);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
