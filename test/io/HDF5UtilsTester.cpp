#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <memory>

#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/io/XFile.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/MPITestUtils.h>

using namespace std;
using namespace xolotl;
using namespace io;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

// Initialize MPI before running any tests; finalize it running all tests.
BOOST_GLOBAL_FIXTURE(MPIFixture);

/**
 * This suite is responsible for testing the HDF5 classes.
 */
BOOST_AUTO_TEST_SUITE(HDF5_testSuite)

/**
 * Create a faux network bound vector.
 */
XFile::NetworkGroup::NetworkBoundsType
createTestNetworkBounds(void)
{
	XFile::NetworkGroup::NetworkBoundsType testBoundsVec{
		{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};
	return testBoundsVec;
}

/**
 * Method checking the writing and reading of the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkIO)
{
	// Determine where we are in the MPI world.
	int commRank = -1;
	int commSize = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	const uint32_t nGridPointsPerRank = 4;
	double const factor = 1.5;

	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 1 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Set the number of grid points and step size
	int nGrid = 5;
	double stepSize = 0.5;
	std::vector<double> grid;
	std::vector<double> temperatures;
	for (int i = 0; i < nGrid + 2; i++) {
		grid.push_back((double)i * stepSize);
		temperatures.push_back(1000.0);
	}

	// Create the network
	using NetworkType = xolotl::core::network::PSIReactionNetwork<
		xolotl::core::network::PSIFullSpeciesList>;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe = opts.getMaxImpurity();
	NetworkType::AmountType maxD = opts.getMaxD();
	NetworkType::AmountType maxT = opts.getMaxT();
	NetworkType network({maxHe, maxD, maxT, maxV, maxI}, grid.size(), opts);
	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);
	// Get the size of the network
	int networkSize = network.getNumClusters();

	// Set the time step number
	int timeStep = 0;

	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;

	// Set the surface information
	XFile::TimestepGroup::Surface1DType iSurface = 3;
	XFile::TimestepGroup::Data1DType nInter = 1.0;
	XFile::TimestepGroup::Data1DType previousFlux = 0.1;
	XFile::TimestepGroup::Data1DType nHe = 1.0;
	XFile::TimestepGroup::Data1DType previousHeFlux = 5.0;
	XFile::TimestepGroup::Data1DType nD = 0.0;
	XFile::TimestepGroup::Data1DType previousDFlux = 0.0;
	XFile::TimestepGroup::Data1DType nT = 0.0;
	XFile::TimestepGroup::Data1DType previousTFlux = 0.0;
	XFile::TimestepGroup::Data1DType nV = 0.5;
	XFile::TimestepGroup::Data1DType previousVFlux = 0.01;

	// Define a faux network composition vector.
	BOOST_TEST_MESSAGE("Creating faux comp vec.");
	auto testBoundsVec = createTestNetworkBounds();

	// Create and populate a Xolotl HDF5 file.
	// We do it in its own scope so that the file is closed
	// when the XFile object goes out of scope.
	const std::string testFileName = "test_basic.h5";
	{
		BOOST_TEST_MESSAGE("Creating file.");
		XFile testFile(testFileName, grid, MPI_COMM_WORLD);
		XFile::NetworkGroup netGroup(testFile, network);
	}

	// Define the concentration dataset size.
	int length = 5;
	int gridPoint = 0;

#if READY
	// Create a vector of concentration for one grid point
	double concArray[length][2];
#endif // READY

	// Define our part of the concentration dataset.
	uint32_t baseX = commRank * nGridPointsPerRank;
	XFile::TimestepGroup::Concs1DType myConcs(nGridPointsPerRank);
	for (auto i = 0; i < nGridPointsPerRank; ++i) {
		auto currNumItems = baseX + i;
		for (auto j = 0; j < currNumItems; ++j) {
			myConcs[i].emplace_back(baseX + j, factor * (baseX + j));
		}
	}

	// Open the test file again to add concentrations.
	// Open it again to add the concentrations.
	// Again, we do this in its own scope so that the file
	// is closed once the object goes out of scope.
	{
		BOOST_TEST_MESSAGE("Opening test file to add a timestep");
		XFile testFile(
			testFileName, MPI_COMM_WORLD, XFile::AccessMode::OpenReadWrite);

		// Add a TimestepGroup.
		auto concGroup = testFile.getGroup<XFile::ConcentrationGroup>();
		BOOST_REQUIRE(concGroup);
		auto tsGroup = concGroup->addTimestepGroup(
			timeStep, currentTime, previousTime, currentTimeStep);

		std::vector<double> nSurf = {nHe, nD, nT, nV};
		std::vector<double> previousSurfFlux = {
			previousHeFlux, previousDFlux, previousTFlux, previousVFlux};
		std::vector<std::string> surfNames = {
			"Helium", "Deuterium", "Tritium", "Vacancy"};
		// Write the surface information
		tsGroup->writeSurface1D(
			iSurface, nInter, previousFlux, nSurf, previousSurfFlux, surfNames);

		std::vector<double> nBulk = {nHe, nV};
		std::vector<double> previousBulkFlux = {previousHeFlux, previousVFlux};
		std::vector<std::string> bulkNames = {"Helium", "Vacancy"};
		// Write the bulk information
		tsGroup->writeBottom1D(nBulk, previousBulkFlux, bulkNames);

#if READY
		// Fill it
		for (int i = 0; i < length; i++) {
			// Fill the concArray
			concArray[i][0] = (double)i;
			concArray[i][1] = (double)i * 10.0 - 5.0;
		}
#endif // READY

		// Write the concentrations in the HDF5 file
#if READY
		tsGroup->writeConcentrationDataset(length, concArray, gridPoint);
#endif // READY
		tsGroup->writeConcentrations(testFile, baseX, myConcs);
	}

	// Now check the test file's contents.
	{
		BOOST_TEST_MESSAGE("Opening test file to check its contents.");
		XFile testFile(
			testFileName, MPI_COMM_WORLD, XFile::AccessMode::OpenReadOnly);

		// Read the header of the written file
		BOOST_TEST_MESSAGE("Checking test file header.");
		int nx = 0, ny = 0, nz = 0;
		double hx = 0.0, hy = 0.0, hz = 0.0;
		auto headerGroup = testFile.getGroup<XFile::HeaderGroup>();
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
		auto concGroup = testFile.getGroup<XFile::ConcentrationGroup>();
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

		// Read the surface information
		BOOST_TEST_MESSAGE(
			"Checking test file last time step surface position.");
		BOOST_REQUIRE_EQUAL(tsGroup->readSurface1D(), iSurface);
		BOOST_REQUIRE_CLOSE(
			tsGroup->readData1D("nInterstitial"), nInter, 0.0001);
		BOOST_REQUIRE_CLOSE(
			tsGroup->readData1D("previousFluxI"), previousFlux, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("nHeliumSurf"), nHe, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("previousFluxHeliumSurf"),
			previousHeFlux, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("nDeuteriumSurf"), nD, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("previousFluxDeuteriumSurf"),
			previousDFlux, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("nTritiumSurf"), nT, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("previousFluxTritiumSurf"),
			previousTFlux, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("nVacancySurf"), nV, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("previousFluxVacancySurf"),
			previousVFlux, 0.0001);

		// Read the bulk information
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("nHeliumBulk"), nHe, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("previousFluxHeliumBulk"),
			previousHeFlux, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("nVacancyBulk"), nV, 0.0001);
		BOOST_REQUIRE_CLOSE(tsGroup->readData1D("previousFluxVacancyBulk"),
			previousVFlux, 0.0001);

		// Read the network of the written file
		BOOST_TEST_MESSAGE("Checking test file last time step network.");
		auto networkGroup = testFile.getGroup<XFile::NetworkGroup>();
		BOOST_REQUIRE(networkGroup);
		int networkSize = networkGroup->readNetworkSize();
		constexpr auto speciesRange = network.getSpeciesRange();
		// Loop on the clusters
		for (int i = 0; i < networkSize; i++) {
			// Get the i-th reactant in the network
			auto cluster = network.getClusterCommon(i);
			// Open the cluster group
			XFile::ClusterGroup clusterGroup(*networkGroup, i);

			// Check the attributes
			double formationEnergy = 0.0, migrationEnergy = 0.0,
				   diffusionFactor = 0.0;
			auto comp = clusterGroup.readCluster(
				formationEnergy, migrationEnergy, diffusionFactor);
			// Check the formation energy
			BOOST_REQUIRE_EQUAL(cluster.getFormationEnergy(), formationEnergy);
			// Check the migration energy
			BOOST_REQUIRE_EQUAL(cluster.getMigrationEnergy(), migrationEnergy);
			// Check the diffusion factor
			BOOST_REQUIRE_EQUAL(cluster.getDiffusionFactor(), diffusionFactor);

			// Check the bounds
			const auto& clReg = network.getCluster(i, plsm::onHost).getRegion();
			NetworkType::Composition lo = clReg.getOrigin();
			NetworkType::Composition hi = clReg.getUpperLimitPoint();
			for (auto j : speciesRange) {
				BOOST_REQUIRE_EQUAL(comp[2 * j()], lo[j]);
				BOOST_REQUIRE_EQUAL(comp[2 * j() + 1], hi[j] - 1);
			}
		}

		// If the HDF5 file contains initial concentrations
		// TODO Doesn't it always contain initial concentrations if
		// it has those earlier items?  Or did I put them on timesteps
		// and they're supposed to be on the concentration group?
		// I suppose it is that one can have a timestep but not
		// have any concentration datasets.
		if (tsGroup) {
			BOOST_TEST_MESSAGE("Testing grid point concentrations.");

#if READY
			// Read the concentrations at the given grid point
			auto returnedVector = tsGroup->readGridPoint(gridPoint);

			// Check the size of the vector
			BOOST_REQUIRE_EQUAL(returnedVector.size(), length);
			// Check the values
			for (unsigned int i = 0; i < returnedVector.size(); i++) {
				BOOST_REQUIRE_CLOSE(
					returnedVector.at(i).at(0), concArray[i][0], 0.0001);
				BOOST_REQUIRE_CLOSE(
					returnedVector.at(i).at(1), concArray[i][1], 0.0001);
			}
#endif // READY

			BOOST_TEST_MESSAGE("Testing grid point concentrations.");

			// Read and check our part of the concentrations.
			auto readConcs = tsGroup->readConcentrations(
				testFile, baseX, nGridPointsPerRank);
			BOOST_REQUIRE_EQUAL(readConcs.size(), myConcs.size());
			for (auto ptIdx = 0; ptIdx < nGridPointsPerRank; ++ptIdx) {
				BOOST_REQUIRE_EQUAL(
					readConcs[ptIdx].size(), myConcs[ptIdx].size());
				for (auto i = 0; i < myConcs[ptIdx].size(); ++i) {
					BOOST_REQUIRE_EQUAL(
						readConcs[ptIdx][i].first, myConcs[ptIdx][i].first);
					BOOST_REQUIRE_CLOSE(readConcs[ptIdx][i].second,
						myConcs[ptIdx][i].second, 0.0001);
				}
			}
		}
	}
}

/**
 * Method checking the writing and reading of the surface position specifically
 * in the case of a 2D grid.
 */
BOOST_AUTO_TEST_CASE(checkSurface2D)
{
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
			grid.push_back((double)i * stepSize);

		XFile testFile(testFileName, grid, MPI_COMM_WORLD);
	}

	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;

	// Define the 2D surface information.
	XFile::TimestepGroup::Surface2DType iSurface = {2, 3, 2, 0, 5};
	XFile::TimestepGroup::Data2DType nInter = {0.0, 0.0, 0.5, 0.6, 0.5};
	XFile::TimestepGroup::Data2DType previousIFlux = {0.0, 0.1, 3.0, -1.0, 5.0};

	// Open the file to add concentrations.
	// Done in its own scope so that it closes when the
	// object goes out of scope.
	{
		BOOST_TEST_MESSAGE("Adding 2D timestep group");

		XFile testFile(
			testFileName, MPI_COMM_WORLD, XFile::AccessMode::OpenReadWrite);

		// Set the time step number
		int timeStep = 0;

		// Add the concentration sub group
		auto concGroup = testFile.getGroup<XFile::ConcentrationGroup>();
		BOOST_REQUIRE(concGroup);
		auto tsGroup = concGroup->addTimestepGroup(
			timeStep, currentTime, previousTime, currentTimeStep);
		BOOST_REQUIRE(tsGroup);

		// Write the surface position
		tsGroup->writeSurface2D(iSurface, nInter, previousIFlux);
	}

	// Read the file to check the values we wrote.
	{
		BOOST_TEST_MESSAGE("Opening 2D file to check its contents.");

		XFile testFile(
			testFileName, MPI_COMM_WORLD, XFile::AccessMode::OpenReadOnly);

		// Access the last written timestep group.
		auto concGroup = testFile.getGroup<XFile::ConcentrationGroup>();
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
		auto previousIFlux = tsGroup->readData2D("previousFluxI");
		// Check all the values
		for (int i = 0; i < previousIFlux.size(); i++) {
			BOOST_REQUIRE_CLOSE(previousIFlux[i], previousIFlux[i], 0.0001);
		}
	}
}

/**
 * Method checking the writing and reading of the surface position specifically
 * in the case of a 3D grid.
 */
BOOST_AUTO_TEST_CASE(checkSurface3D)
{
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
			grid.push_back((double)i * stepSize);

		XFile testFile(testFileName, grid, MPI_COMM_WORLD);
	}

	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;

	// Set the 3D surface information.
	XFile::TimestepGroup::Surface3DType iSurface = {
		{2, 4, 1, 0, 5}, {2, 3, 2, 0, 5}, {6, 1, 2, 3, 2}};
	XFile::TimestepGroup::Data3DType nInter = {{0.0, 0.0, 0.0, 0.0, 0.0},
		{2.0, 3.0, 2.0, 0.0, 0.5}, {0.0, 0.0, 0.0, 0.0, 0.0}};
	XFile::TimestepGroup::Data3DType previousIFlux = {{0.0, 0.0, 0.0, 0.0, 0.0},
		{-2.0, 3.0, 2.0, 0.0, -0.5}, {0.0, 0.0, 0.0, 0.0, 0.0}};

	// Open test file to add timestep group with concentrations.
	// Done in its own scope so that it closes when the
	// object goes out of scope.
	{
		BOOST_TEST_MESSAGE("Adding 3D timestep group");

		XFile testFile(
			testFileName, MPI_COMM_WORLD, XFile::AccessMode::OpenReadWrite);

		// Set the time step number
		int timeStep = 0;

		// Add the concentration sub group
		auto concGroup = testFile.getGroup<XFile::ConcentrationGroup>();
		BOOST_REQUIRE(concGroup);

		auto tsGroup = concGroup->addTimestepGroup(
			timeStep, currentTime, previousTime, currentTimeStep);
		BOOST_REQUIRE(tsGroup);

		// Write the surface position
		tsGroup->writeSurface3D(iSurface, nInter, previousIFlux);
	}

	// Check contents of file we wrote.
	{
		BOOST_TEST_MESSAGE("Opening 3D file to check its contents.");
		XFile testFile(
			testFileName, MPI_COMM_WORLD, XFile::AccessMode::OpenReadOnly);

		// Access the last written timestep group.
		auto concGroup = testFile.getGroup<XFile::ConcentrationGroup>();
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
		auto previousIFlux = tsGroup->readData3D("previousFluxI");
		// Check all the values
		for (int i = 0; i < previousIFlux.size(); i++) {
			for (int j = 0; j < previousIFlux[0].size(); j++) {
				BOOST_REQUIRE_CLOSE(
					previousIFlux[i][j], previousIFlux[i][j], 0.0001);
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
