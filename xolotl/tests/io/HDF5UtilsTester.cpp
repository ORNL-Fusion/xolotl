#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <boost/test/framework.hpp>
#include <PSIClusterReactionNetwork.h>
#include <DummyHandlerRegistry.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>
#include <Options.h>
#include "xolotlCore/io/XFile.h"

using namespace std;
using namespace xolotlCore;

struct MPIFixture {
    MPIFixture(void) {
        BOOST_TEST_MESSAGE("Initializing MPI");
        auto& mts = boost::unit_test::framework::master_test_suite();
        MPI_Init(&mts.argc, &mts.argv);
        BOOST_TEST_MESSAGE("Done initializing MPI");
    }

    ~MPIFixture(void) {
        BOOST_TEST_MESSAGE("Finalizing MPI");
        MPI_Finalize();
        BOOST_TEST_MESSAGE("Done finalizing MPI");
    }
};
BOOST_GLOBAL_FIXTURE(MPIFixture);

/**
 * This suite is responsible for testing the HDF5 classes.
 */
BOOST_AUTO_TEST_SUITE(HDF5_testSuite)

/**
 * Create a faux network composition vector.
 */
xolotlCore::XFile::HeaderGroup::NetworkCompsType createTestNetworkComps(void) {

    std::vector<std::vector<int>> testCompsVec {
        { 1, 2, 3, 4 },
        { 5, 6, 7, 8 },
        { 9, 10, 11, 12 }
    };
    return testCompsVec;
}


/**
 * Method checking the writing and reading of the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkIO) {
    BOOST_TEST_MESSAGE("Creating network loader.");
	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

    BOOST_TEST_MESSAGE("Creating network load options.");
	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Set the time step number
	int timeStep = 0;

	// Set the number of grid points and step size
	int nGrid = 5;
	double stepSize = 0.5;
	std::vector<double> grid;
	for (int i = 0; i < nGrid + 2; i++)
		grid.push_back((double) i * stepSize);
	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;

	// Set the surface information
    xolotlCore::XFile::TimestepGroup::Surface1DType iSurface = 3;
    xolotlCore::XFile::TimestepGroup::Data1DType nInter = 1.0;
    xolotlCore::XFile::TimestepGroup::Data1DType previousFlux = 0.1;

    // Define a faux network composition vector.
    BOOST_TEST_MESSAGE("Creating faux comp vec.");
    auto testCompsVec = createTestNetworkComps();

    // Create and populate a Xolotl HDF5 file.
    // We do it in its own scope so that the file is closed
    // when the XFile object goes out of scope.
    const std::string testFileName = "test_basic.h5";
    {
        BOOST_TEST_MESSAGE("Creating file.");
        xolotlCore::XFile testFile(testFileName,
                                    grid,
                                    testCompsVec, 
                                    filename);
        std::cout << "Done creating file" << std::endl;
    }


    // Define the concentration dataset size.
    int length = 5;
    int gridPoint = 0;

    // Create a vector of concentration for one grid point
    double concArray[length][2];

    // Open the test file again to add concentrations.
	// Open it again to add the concentrations.
    // Again, we do this in its own scope so that the file
    // is closed once the object goes out of scope.
    {
        std::cout << "Opening file" << std::endl;
        xolotlCore::XFile testFile(testFileName,
                xolotlCore::XFile::AccessMode::OpenReadWrite);

        // Add a TimestepGroup.
        auto concGroup = testFile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
        BOOST_REQUIRE(concGroup);
        auto tsGroup = concGroup->addTimestepGroup(timeStep, currentTime,
                                    previousTime, currentTimeStep);

        // Write the surface position
        tsGroup->writeSurface1D(iSurface, nInter, previousFlux);

        // Fill it
        for (int i = 0; i < length; i++) {
            // Fill the concArray
            concArray[i][0] = (double) i;
            concArray[i][1] = (double) i * 10.0 - 5.0;
        }

        // Write the concentrations in the HDF5 file
        tsGroup->writeConcentrationDataset(length, concArray, gridPoint);
        std::cout << "Done writing concentrations" << std::endl;
    }

    // Now check the test file's contents.
    {
        BOOST_TEST_MESSAGE("Opening test file to check its contents.");
        xolotlCore::XFile testFile(testFileName,
                xolotlCore::XFile::AccessMode::OpenReadOnly);

        // Read the header of the written file
        BOOST_TEST_MESSAGE("Checking test file header.");
        int nx = 0, ny = 0, nz = 0;
        double hx = 0.0, hy = 0.0, hz = 0.0;
        auto headerGroup = testFile.getGroup<xolotlCore::XFile::HeaderGroup>();
        BOOST_REQUIRE(headerGroup);
        headerGroup->read(nx, hx, ny, hy, nz, hz);

        // Check the obtained values
        BOOST_REQUIRE_EQUAL(nx, nGrid);
        BOOST_REQUIRE_CLOSE(hx, stepSize, 0.0001);
        BOOST_REQUIRE_EQUAL(ny, 0);
        BOOST_REQUIRE_CLOSE(hy, 0.0, 0.0001);
        BOOST_REQUIRE_EQUAL(nz, 0);
        BOOST_REQUIRE_CLOSE(hz, 0.0, 0.0001);

        // Access the last-written timestep group.
        BOOST_TEST_MESSAGE("Checking test file last time step.");
        auto concGroup = testFile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
        BOOST_REQUIRE(concGroup);
        auto tsGroup = concGroup->getLastTimestepGroup();
        BOOST_REQUIRE(tsGroup);

        // Read the times
        BOOST_TEST_MESSAGE("Checking test file last time step times.");
        double t = 0.0, dt = 0.0;
        std::tie(t, dt) = tsGroup->readTimes();
        BOOST_REQUIRE_CLOSE(t, currentTime, 0.0001);
        BOOST_REQUIRE_CLOSE(dt, currentTimeStep, 0.0001);
        double previousReadTime = tsGroup->readPreviousTime();
        BOOST_REQUIRE_CLOSE(previousReadTime, previousTime, 0.0001);

        // Read the surface position
        BOOST_TEST_MESSAGE("Checking test file last time step surface position.");
        auto surfacePos = tsGroup->readSurface1D();
        BOOST_REQUIRE_EQUAL(surfacePos, iSurface);
        auto nInterstitial = tsGroup->readData1D("nInterstitial");
        BOOST_REQUIRE_CLOSE(nInterstitial, nInter, 0.0001);
        auto previousIFlux = tsGroup->readData1D("previousIFlux");
        BOOST_REQUIRE_CLOSE(previousIFlux, previousFlux, 0.0001);

        // Read the network of the written file
        BOOST_TEST_MESSAGE("Checking test file last time step network.");
        auto networkGroup = testFile.getGroup<xolotlCore::XFile::NetworkGroup>();
        BOOST_REQUIRE(networkGroup);
        auto networkVector = networkGroup->readNetwork();
        // Get all the reactants
        auto const& reactants = network->getAll();
        // Check the network vector
        for (IReactant& it : reactants) {
            // Get the i-th reactant in the network
            auto& reactant = (PSICluster&) it;
            int id = reactant.getId() - 1;
            // Get the corresponding line from the HDF5 file
            auto line = networkVector.at(id);

            // Check the composition
            auto& composition = reactant.getComposition();
            BOOST_REQUIRE_EQUAL((int ) line[0],
                    composition[toCompIdx(Species::He)]);
            BOOST_REQUIRE_EQUAL((int ) line[1], composition[toCompIdx(Species::D)]);
            BOOST_REQUIRE_EQUAL((int ) line[2], composition[toCompIdx(Species::T)]);
            BOOST_REQUIRE_EQUAL((int ) line[3], composition[toCompIdx(Species::V)]);
            BOOST_REQUIRE_EQUAL((int ) line[4], composition[toCompIdx(Species::I)]);

            // Check the formation energy
            auto formationEnergy = reactant.getFormationEnergy();
            BOOST_REQUIRE_EQUAL(line[5], formationEnergy);

            // Check the migration energy
            double migrationEnergy = reactant.getMigrationEnergy();
            BOOST_REQUIRE_EQUAL(line[6], migrationEnergy);

            // Check the diffusion factor
            double diffusionFactor = reactant.getDiffusionFactor();
            BOOST_REQUIRE_EQUAL(line[7], diffusionFactor);
        }

        // Test the composition vector.
        BOOST_TEST_MESSAGE("Checking test file last time step composition vectors.");
        auto readCompsVec = headerGroup->readNetworkComps();
        BOOST_REQUIRE_EQUAL(readCompsVec.size(), testCompsVec.size());
        for(auto i = 0; i < testCompsVec.size(); ++i) {
            BOOST_REQUIRE_EQUAL(readCompsVec[i].size(), testCompsVec[i].size());
            for(auto j = 0; j < testCompsVec[i].size(); ++j) {
                BOOST_REQUIRE_EQUAL(readCompsVec[i][j], testCompsVec[i][j]);
            }
        }

        // If the HDF5 file contains initial concentrations
        // TODO Doesn't it always contain initial concentrations if
        // it has those earlier items?  Or did I put them on timesteps
        // and they're supposed to be on the concentration group?
        // I suppose it is that one can have a timestep but not
        // have any concentration datasets.
        if(tsGroup) {

            BOOST_TEST_MESSAGE("Testing grid point concentrations.");

            // Read the concentrations at the given grid point
            auto returnedVector = tsGroup->readGridPoint(gridPoint);

            // Check the size of the vector
            BOOST_REQUIRE_EQUAL(returnedVector.size(), length);
            // Check the values
            for (unsigned int i = 0; i < returnedVector.size(); i++) {
                BOOST_REQUIRE_CLOSE(returnedVector.at(i).at(0), concArray[i][0],
                        0.0001);
                BOOST_REQUIRE_CLOSE(returnedVector.at(i).at(1), concArray[i][1],
                        0.0001);
            }
        }
    }
}

/**
 * Method checking the writing and reading of the surface position specifically
 * in the case of a 2D grid.
 */
BOOST_AUTO_TEST_CASE(checkSurface2D) {


    // Create the test HDF5 file.
    // Done in its own scope so that it closes when the 
    // object goes out of scope.
    const std::string testFileName = "test_surface2D.h5";
    {
        BOOST_TEST_MESSAGE("Creating 2D test file");

        // Set the number of grid points and step size
        int nGrid = 5;
        double stepSize = 0.5;
        std::vector<double> grid;
        for (int i = 0; i < nGrid + 2; i++)
            grid.push_back((double) i * stepSize);

        xolotlCore::XFile testFile(testFileName,
                                    grid,
                                    createTestNetworkComps(), 
                                    "");        // no file to copy network from
    }

    // Set the time information
    double currentTime = 0.0001;
    double previousTime = 0.00001;
    double currentTimeStep = 0.000001;

    // Define the 2D surface information.
    xolotlCore::XFile::TimestepGroup::Surface2DType iSurface = { 2, 3, 2, 0, 5 };
    xolotlCore::XFile::TimestepGroup::Data2DType nInter = { 0.0, 0.0, 0.5, 0.6, 0.5 };
    xolotlCore::XFile::TimestepGroup::Data2DType previousFlux = { 0.0, 0.1, 3.0, -1.0, 5.0 };

    // Open the file to add concentrations.
    // Done in its own scope so that it closes when the
    // object goes out of scope.
    {
        BOOST_TEST_MESSAGE("Adding 2D timestep group");

        xolotlCore::XFile testFile(testFileName,
                xolotlCore::XFile::AccessMode::OpenReadWrite);

        // Set the time step number
        int timeStep = 0;

        // Add the concentration sub group
        auto concGroup = testFile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
        BOOST_REQUIRE(concGroup);
        auto tsGroup = concGroup->addTimestepGroup(timeStep, currentTime,
                                            previousTime, currentTimeStep);
        BOOST_REQUIRE(tsGroup);

        // Write the surface position
        tsGroup->writeSurface2D(iSurface, nInter, previousFlux);
    }

    // Read the file to check the values we wrote.
    {
        BOOST_TEST_MESSAGE("Opening 2D file to check its contents.");

        xolotlCore::XFile testFile(testFileName,
                xolotlCore::XFile::AccessMode::OpenReadOnly);

        // Access the last written timestep group.
        auto concGroup = testFile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
        BOOST_REQUIRE(concGroup);
        auto tsGroup = concGroup->getLastTimestepGroup();
        BOOST_REQUIRE(tsGroup);

        // Read the surface position
        auto surfacePos = tsGroup->readSurface2D();
        // Check all the values
        for (int i = 0; i < surfacePos.size(); i++) {
            BOOST_REQUIRE_EQUAL(surfacePos[i], iSurface[i]);
        }

        // Read the interstitial quantity
        auto nInterstitial = tsGroup->readData2D("nInterstitial");
        // Check all the values
        for (int i = 0; i < nInterstitial.size(); i++) {
            BOOST_REQUIRE_CLOSE(nInterstitial[i], nInter[i], 0.0001);
        }

        // Read the interstitial flux
        auto previousIFlux = tsGroup->readData2D("previousIFlux");
        // Check all the values
        for (int i = 0; i < previousIFlux.size(); i++) {
            BOOST_REQUIRE_CLOSE(previousIFlux[i], previousFlux[i], 0.0001);
        }
    }
}

/**
 * Method checking the writing and reading of the surface position specifically
 * in the case of a 3D grid.
 */
BOOST_AUTO_TEST_CASE(checkSurface3D) {

    const std::string testFileName = "test_surface3D.h5";

    // Create test file.
    // Done in its own scope so that it closes when the
    // object goes out of scope.
    {
        BOOST_TEST_MESSAGE("Creating 3D test file.");

        // Set the number of grid points and step size
        int nGrid = 5;
        double stepSize = 0.5;
        std::vector<double> grid;
        for (int i = 0; i < nGrid + 2; i++)
            grid.push_back((double) i * stepSize);

        xolotlCore::XFile testFile(testFileName,
                                    grid,
                                    createTestNetworkComps(),
                                    "");    // no network file to copy from
    }

    // Set the time information
    double currentTime = 0.0001;
    double previousTime = 0.00001;
    double currentTimeStep = 0.000001;

	// Set the 3D surface information.
    xolotlCore::XFile::TimestepGroup::Surface3DType iSurface = {
        { 2, 4, 1, 0, 5 },
        { 2, 3, 2, 0, 5 },
        { 6, 1, 2, 3, 2 }
    };
    xolotlCore::XFile::TimestepGroup::Data3DType nInter = {
        { 0.0, 0.0, 0.0, 0.0, 0.0 },
        {2.0, 3.0, 2.0, 0.0, 0.5 },
        { 0.0, 0.0, 0.0, 0.0, 0.0 }
    };
    xolotlCore::XFile::TimestepGroup::Data3DType previousFlux = {
        { 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -2.0, 3.0, 2.0, 0.0, -0.5 },
        { 0.0, 0.0, 0.0, 0.0, 0.0 }
    };

    // Open test file to add timestep group with concentrations.
    // Done in its own scope so that it closes when the
    // object goes out of scope.
    {
        BOOST_TEST_MESSAGE("Adding 3D timestep group");

        xolotlCore::XFile testFile(testFileName,
                xolotlCore::XFile::AccessMode::OpenReadWrite);

        // Set the time step number
        int timeStep = 0;

        // Add the concentration sub group
        auto concGroup = testFile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
        BOOST_REQUIRE(concGroup);

        auto tsGroup = concGroup->addTimestepGroup(timeStep, currentTime,
                                            previousTime, currentTimeStep);
        BOOST_REQUIRE(tsGroup);

        // Write the surface position
        tsGroup->writeSurface3D(iSurface, nInter, previousFlux);
    }

    // Check contents of file we wrote.
    {
        BOOST_TEST_MESSAGE("Opening 3D file to check its contents.");
        xolotlCore::XFile testFile(testFileName,
                xolotlCore::XFile::AccessMode::OpenReadOnly);

        // Access the last written timestep group.
        auto concGroup = testFile.getGroup<xolotlCore::XFile::ConcentrationGroup>();
        BOOST_REQUIRE(concGroup);
        auto tsGroup = concGroup->getLastTimestepGroup();
        BOOST_REQUIRE(tsGroup);

        // Read the surface position
        auto surfacePos = tsGroup->readSurface3D();
        // Check all the values
        for (int i = 0; i < surfacePos.size(); i++) {
            for (int j = 0; j < surfacePos[0].size(); j++) {
                BOOST_REQUIRE_EQUAL(surfacePos[i][j], iSurface[i][j]);
            }
        }
        auto nInterstitial = tsGroup->readData3D("nInterstitial");
        // Check all the values
        for (int i = 0; i < nInterstitial.size(); i++) {
            for (int j = 0; j < nInterstitial[0].size(); j++) {
                BOOST_REQUIRE_CLOSE(nInterstitial[i][j], nInter[i][j], 0.0001);
            }
        }
        auto previousIFlux = tsGroup->readData3D("previousIFlux");
        // Check all the values
        for (int i = 0; i < previousIFlux.size(); i++) {
            for (int j = 0; j < previousIFlux[0].size(); j++) {
                BOOST_REQUIRE_CLOSE(previousIFlux[i][j], previousFlux[i][j],
                        0.0001);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
