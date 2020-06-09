#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <xolotl/core/network/NEReactionNetwork.h>

using namespace std;
using namespace xolotl::core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the NE network.
 */
BOOST_AUTO_TEST_SUITE(NENetwork_testSuite)

BOOST_AUTO_TEST_CASE(fullyRefined) {

	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=20 0 0 0 0" << std::endl << "process=reaction"
			<< std::endl;
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
	opts.readParams(argc, argv);

	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;
	NetworkType network( { opts.getMaxImpurity() }, 1, opts);

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 20);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 20);
	// TODO: check it is within a given range?
	BOOST_REQUIRE(network.getDeviceMemorySize() > 0);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.574, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.047279806, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 20);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 1);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 1);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 1);

	// Check the single vacancy
	auto vacancy = network.getSingleVacancy();
	BOOST_REQUIRE(vacancy.getId() == network.invalidIndex());

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
			17, 18, 19 };
	knownDFill[1] = { 1, 0, 2 };
	knownDFill[2] = { 2, 0, 1, 3 };
	knownDFill[3] = { 3, 0, 2, 4 };
	knownDFill[4] = { 4, 0, 3, 5 };
	knownDFill[5] = { 5, 0, 4, 6 };
	knownDFill[6] = { 6, 0, 5, 7 };
	knownDFill[7] = { 7, 0, 6, 8 };
	knownDFill[8] = { 8, 0, 7, 9 };
	knownDFill[9] = { 9, 0, 8, 10 };
	knownDFill[10] = { 10, 0, 9, 11 };
	knownDFill[11] = { 11, 0, 10, 12 };
	knownDFill[12] = { 12, 0, 11, 13 };
	knownDFill[13] = { 13, 0, 12, 14 };
	knownDFill[14] = { 14, 0, 13, 15 };
	knownDFill[15] = { 15, 0, 14, 16 };
	knownDFill[16] = { 16, 0, 15, 17 };
	knownDFill[17] = { 17, 0, 16, 18 };
	knownDFill[18] = { 18, 0, 17, 19 };
	knownDFill[19] = { 19, 0, 18 };
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 94);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto &row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				BOOST_REQUIRE_EQUAL(row[j], knownDFill[i][j]);
			}
		}
	}

	// Set temperatures
	std::vector<double> temperatures = { 1000.0 };
	network.setTemperatures(temperatures);
	network.syncClusterDataOnHost();
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 5.4736558e-6, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
	Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(network.getLeftSideRate(dConcs, 0, gridId), 7.98584e-5,
			0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
	Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = { -8.5332e-05, 2.45933e-06, -2.38186e-07,
			-1.8962e-07, -1.60127e-07, -1.39969e-07, -1.25163e-07, -1.13743e-07,
			-1.04618e-07, -9.71289e-08, -9.08512e-08, -8.54991e-08,
			-8.08717e-08, -7.68238e-08, -7.32473e-08, -7.006e-08, -6.71984e-08,
			-6.46122e-08, -6.22614e-08, 4.8543e-06, 0 };
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = { -9.62794e-05, -3.01432e-06,
			-3.25251e-06, -3.44213e-06, -3.60226e-06, -3.74222e-06,
			-3.86739e-06, -3.98113e-06, -4.08575e-06, -4.18288e-06,
			-4.27373e-06, -4.35923e-06, -4.4401e-06, -4.51692e-06, -4.59017e-06,
			-4.66023e-06, -4.72743e-06, -4.79204e-06, -4.8543e-06, 3.50475e-32,
			-3.01432e-06, 7.93299e-06, 5.30831e-15, -3.25251e-06, -2.38186e-07,
			3.01432e-06, 3.14804e-16, -3.44213e-06, -1.8962e-07, 3.25251e-06,
			5.84363e-17, -3.60226e-06, -1.60127e-07, 3.44213e-06, 3.42327e-17,
			-3.74222e-06, -1.39969e-07, 3.60226e-06, 1.09411e-18, -3.86739e-06,
			-1.25163e-07, 3.74222e-06, 3.54301e-19, -3.98113e-06, -1.13743e-07,
			3.86739e-06, 3.45219e-22, -4.08575e-06, -1.04618e-07, 3.98113e-06,
			2.05723e-20, -4.18288e-06, -9.71289e-08, 4.08575e-06, 3.69419e-21,
			-4.27373e-06, -9.08512e-08, 4.18288e-06, 1.14015e-23, -4.35923e-06,
			-8.54991e-08, 4.27373e-06, 6.39174e-25, -4.4401e-06, -8.08717e-08,
			4.35923e-06, 2.03997e-25, -4.51692e-06, -7.68238e-08, 4.4401e-06,
			3.64005e-26, -4.59017e-06, -7.32473e-08, 4.51692e-06, 2.03305e-27,
			-4.66023e-06, -7.006e-08, 4.59017e-06, 6.4677e-28, -4.72743e-06,
			-6.71984e-08, 4.66023e-06, 6.32493e-30, -4.79204e-06, -6.46122e-08,
			4.72743e-06, 3.52376e-31, -4.8543e-06, -6.22614e-08, 4.79204e-06,
			3.50475e-32, -3.50475e-32, 4.8543e-06, 4.8543e-06 };
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	for (NetworkType::IndexType i = 0; i < nPartials; i++) {
		BOOST_REQUIRE_CLOSE(hPartials[i], knownPartials[i], 0.01);
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::Xe] = 1;
	auto cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 3.6298e-7,
			0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.0, 0.01);
	const auto &reg = cluster.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 2);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Xe] = 10;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 9);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.61702, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 46.9, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(cluster.getMigrationEnergy(),
			std::numeric_limits<double>::infinity());
	const auto &regBis = cluster.getRegion();
	lo = regBis.getOrigin();
	hi = regBis.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 10);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 11);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
}

BOOST_AUTO_TEST_SUITE_END()
