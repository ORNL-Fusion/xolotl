#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <PSIClusterReactionNetwork.h>
#include <DummyHandlerRegistry.h>
#include <PSIClusterNetworkLoader.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HDF5Utils
 */
BOOST_AUTO_TEST_SUITE(HDF5Utils_testSuite)

/**
 * Method checking the writing and reading of the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkOI) {

	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	PSIClusterNetworkLoader loader =
			PSIClusterNetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten.txt");
	string filename = sourceDir + pathToFile;
	// Create the network stream
	shared_ptr<istream> networkStream;
	networkStream = make_shared<ifstream>(filename);
	// Read the buffer of the stream
	auto bufferSS = make_shared<stringstream>();
	(*bufferSS) << networkStream->rdbuf();
	// Give the network stream to the network loader
	loader.setInputstream(bufferSS);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Get the size of the network
	int networkSize = network->size();
	// Set the time step number
	int timeStep = 0;
	// Initialize the HDF5 file
	HDF5Utils::initializeFile("test.h5", networkSize);

	// Set the number of grid points and step size
	int nGrid = 5;
	double stepSize = 0.5;
	// Set the time information
	double currentTime = 0.0001;
	double currentTimeStep = 0.000001;
	// Write the header in the HDF5 file
	HDF5Utils::fillHeader(nGrid, stepSize);

	// Write the network in the HDF5 file
	HDF5Utils::fillNetwork(network);

	// Finalize the HDF5 file
	HDF5Utils::finalizeFile();

	// Open it again to add the concentrations
	HDF5Utils::openFile("test.h5");

	// Add the concentration sub group
	HDF5Utils::addConcentrationSubGroup(timeStep, networkSize, currentTime, currentTimeStep);

	// Add the concentration dataset
	int length = 5;
	int gridPoint = 0;
	HDF5Utils::addConcentrationDataset(length, gridPoint);

	// Create a vector of concentration for one grid point
	std::vector< std::vector<double> > concVector;
	// Fill it
	for (int i = 0; i < length; i++) {
		// Create the concentration vector for this cluster
		std::vector<double> conc;
		conc.push_back((double) i);
		conc.push_back((double) i * 10.0 - 5.0);

		// Add it to the main vector
		concVector.push_back(conc);
	}

	// Write the concentrations in the HDF5 file
	HDF5Utils::fillConcentrations(concVector, gridPoint);

	// Close the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	// Read the header of the written file
	int nx = 0, ny = 0, nz = 0;
	double hx = 0.0, hy = 0.0, hz = 0.0;
	HDF5Utils::readHeader("test.h5", nx, hx, ny, hy, nz, hz);
	// Check the obtained values
	BOOST_REQUIRE_EQUAL(nx, nGrid);
	BOOST_REQUIRE_CLOSE(hx, stepSize, 0.0001);
	BOOST_REQUIRE_EQUAL(ny, 0);
	BOOST_REQUIRE_CLOSE(hy, 0.0, 0.0001);
	BOOST_REQUIRE_EQUAL(nz, 0);
	BOOST_REQUIRE_CLOSE(hz, 0.0, 0.0001);

	// Read the times
	double t = 0.0, dt = 0.0;
	HDF5Utils::readTimes("test.h5", 0, t, dt);
	BOOST_REQUIRE_EQUAL(t, currentTime);
	BOOST_REQUIRE_EQUAL(dt, currentTimeStep);

	// Read the network of the written file
	auto networkVector = HDF5Utils::readNetwork("test.h5");
	// Get all the reactants
	auto reactants = network->getAll();
	// Check the network vector
	for (int i = 0; i < networkSize; i++) {
		// Get the i-th reactant in the network
		auto reactant = (PSICluster *) reactants->at(i);
		int id = reactant->getId() - 1;
		// Get the corresponding line from the HDF5 file
		auto line = networkVector.at(id);

		// Check the composition
		auto composition = reactant->getComposition();
		BOOST_REQUIRE_EQUAL((int) line[0], composition["He"]);
		BOOST_REQUIRE_EQUAL((int) line[1], composition["V"]);
		BOOST_REQUIRE_EQUAL((int) line[2], composition["I"]);

		// Check the formation energy
		auto formationEnergy = reactant->getFormationEnergy();
		BOOST_REQUIRE_EQUAL(line[3], formationEnergy);

		// Check the migration energy
		double migrationEnergy = reactant->getMigrationEnergy();
		BOOST_REQUIRE_EQUAL(line[4], migrationEnergy);

		// Check the diffusion factor
		double diffusionFactor = reactant->getDiffusionFactor();
		BOOST_REQUIRE_EQUAL(line[5], diffusionFactor);
	}

	// If the HDF5 file contains initial concentrations
	int tempTimeStep = -2;
	if (HDF5Utils::hasConcentrationGroup("test.h5", tempTimeStep)) {
		// Read the concentrations at the given grid point
		double newConcentrations[networkSize];
		double * newConc = &newConcentrations[0];
		auto returnedVector = HDF5Utils::readGridPoint("test.h5", tempTimeStep, gridPoint);

		// Check the size of the vector
		BOOST_REQUIRE_EQUAL(returnedVector.size(), concVector.size());
		// Check the values
		for (int i = 0; i < returnedVector.size(); i++) {
			BOOST_REQUIRE_CLOSE(returnedVector.at(i).at(0), concVector.at(i).at(0), 0.0001);
			BOOST_REQUIRE_CLOSE(returnedVector.at(i).at(1), concVector.at(i).at(1), 0.0001);
		}
	}

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
