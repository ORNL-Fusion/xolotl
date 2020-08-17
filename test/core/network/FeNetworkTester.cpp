#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/FeReactionNetwork.h>

using namespace std;
using namespace xolotl::core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the Fe network.
 */
BOOST_AUTO_TEST_SUITE(FeNetwork_testSuite)

BOOST_AUTO_TEST_CASE(fullyRefined)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=3 0 0 3 1" << std::endl
			  << "process=reaction" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char** argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	opts.readParams(argc, argv);

	using NetworkType = FeReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;
	NetworkType network({(NetworkType::AmountType)opts.getMaxImpurity(),
							(NetworkType::AmountType)opts.getMaxV(),
							(NetworkType::AmountType)opts.getMaxI()},
		1, opts);

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 16);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 16);
	// TODO: check it is within a given range?
	auto deviceMemorySize = network.getDeviceMemorySize();
	BOOST_REQUIRE(deviceMemorySize > 90000);
	BOOST_REQUIRE(deviceMemorySize < 110000);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.287, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.01181995, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 16);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 2);

	// Check the single vacancy
	auto vacancy = network.getSingleVacancy();
	BOOST_REQUIRE_EQUAL(vacancy.getId(), 1);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 4, 8, 12};
	knownDFill[1] = {1, 0, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 3, 7, 11, 15};
	knownDFill[2] = {2, 0, 3, 1, 4, 5, 8, 9, 12, 13, 6};
	knownDFill[3] = {3, 0, 1, 2, 4, 8, 12, 7};
	knownDFill[4] = {4, 0, 5, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	knownDFill[5] = {5, 0, 6, 1, 4, 2, 8, 9};
	knownDFill[6] = {6, 0, 7, 1, 5, 2, 4, 8, 10};
	knownDFill[7] = {7, 0, 1, 6, 2, 5, 3, 4, 8, 11};
	knownDFill[8] = {8, 0, 9, 1, 2, 3, 4, 5, 6, 7, 12};
	knownDFill[9] = {9, 0, 10, 1, 8, 2, 4, 5, 13};
	knownDFill[10] = {10, 0, 11, 1, 9, 2, 8, 4, 6, 14};
	knownDFill[11] = {11, 0, 1, 10, 2, 9, 3, 8, 4, 7, 15};
	knownDFill[12] = {12, 0, 13, 1, 2, 3, 4, 8};
	knownDFill[13] = {13, 0, 14, 1, 12, 2, 4, 9, 5, 8};
	knownDFill[14] = {14, 0, 15, 1, 13, 2, 12, 4, 10, 6, 8};
	knownDFill[15] = {15, 0, 1, 14, 2, 13, 3, 12, 4, 11, 7, 8};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 176);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	network.setTemperatures(temperatures);
	network.syncClusterDataOnHost();
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 751630489466.0, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 1), 12.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 3), 12.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::I, 1), 0.14131, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::He, 1), 18.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 94327839814, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-9.43288e+10, -5.16275e+11, -5.63061e+11,
		-6.10778e+11, -4.75639e+12, -1.45044e+11, -1.55625e+11, -1.69757e+11,
		-7.60177e+11, 1.45774e+11, 1.57282e+11, 1.59414e+11, 2.6994e+11,
		5.20973e+11, 5.63078e+11, 5.86753e+11, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-9.43288e+10, -7.01816e+09,
		-7.91224e+09, -1.03421e+10, -6.86896e+09, -7.76165e+09, -8.38785e+09,
		-6.86896e+09, -7.76165e+09, -8.38785e+09, -6.86896e+09, -7.76165e+09,
		-8.38785e+09, 1.10364e-10, 0.00140902, 38.3638, -5.30342e+11,
		8.94081e+08, 1.2078e+10, -2.76652e+11, 7.84608e+08, -1.68579e+08,
		-1.45181e+11, -1.49194e+08, -1.68588e+08, -9.90247e+10, -1.49198e+08,
		-1.68588e+08, 1.39186e+08, 17448.6, 26.9324, 0.13587, -5.73858e+11,
		2.42983e+09, 1.04813e+10, 2.77618e+08, -2.9965e+11, -1.50587e+08,
		-1.56675e+11, -1.50587e+08, -1.06607e+11, -1.50587e+08, 1.75351e+07,
		-6.1113e+11, -1.03421e+10, 3.19175e+08, 3.19175e+08, -3.18416e+11,
		-1.67483e+11, -1.1475e+11, 3.29555e+07, -6.26796e+12, 6.86896e+09,
		-2.68616e+11, -2.76652e+11, -2.9965e+11, -3.18416e+11, -2.99407e+11,
		-3.1553e+11, -5.83904e+11, -2.76419e+11, -2.99425e+11, -3.15563e+11,
		4.50762e+08, 1.47931e-12, 1.60243e-12, 1.6888e-12, -4.29458e+11,
		8.92693e+08, 7.76166e+09, 2.76503e+11, 2.32972e+08, -1.50587e+08,
		-1.44936e+11, 1.47931e-12, -4.63812e+11, 6.26203e+08, 8.38787e+09,
		-1.93899e+07, 1.49198e+08, 2.9965e+11, 2.25415e+08, -1.56439e+11,
		1.60243e-12, -4.88492e+11, -8.38785e+09, 1.68588e+08, 1.68588e+08,
		1.50587e+08, 1.50587e+08, 3.18416e+11, 2.85278e+09, -1.64509e+11,
		1.6888e-12, -1.51913e+12, 6.86896e+09, 6.86896e+09, -1.45181e+11,
		-1.56675e+11, -1.67483e+11, 9.19357e+11, -1.44936e+11, -1.56439e+11,
		-1.64509e+11, 4.50762e+08, -2.83588e+11, 8.92693e+08, 7.76165e+09,
		1.45032e+11, 1.45181e+11, -1.50587e+08, 0, 2.76419e+11, 1.47931e-12,
		-3.07355e+11, 6.26203e+08, 8.38785e+09, -1.93899e+07, 1.49198e+08,
		1.56675e+11, 1.56675e+11, 0, 2.99425e+11, 1.60243e-12, -3.23951e+11,
		-8.38785e+09, 1.68588e+08, 1.68588e+08, 1.50587e+08, 1.50587e+08,
		1.67483e+11, 1.67483e+11, 0, 3.15563e+11, 1.6888e-12, -3.20832e+11,
		6.86896e+09, 6.86896e+09, -9.90247e+10, -1.06607e+11, -1.1475e+11,
		5.83904e+11, 5.83904e+11, -7.16874e+09, 8.92693e+08, 7.76165e+09,
		9.88755e+10, 9.90247e+10, -1.50587e+08, 2.76419e+11, 2.76419e+11,
		1.44936e+11, 1.44936e+11, -7.93024e+09, 6.26203e+08, 8.38785e+09,
		-1.93899e+07, 1.49198e+08, 1.06607e+11, 1.06607e+11, 2.99425e+11,
		2.99425e+11, 1.56439e+11, 1.56439e+11, -8.38785e+09, -8.38785e+09,
		1.68588e+08, 1.68588e+08, 1.50587e+08, 1.50587e+08, 1.1475e+11,
		1.1475e+11, 3.15563e+11, 3.15563e+11, 1.64509e+11, 1.64509e+11};
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					knownPartials[startingIdx + j], 0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::He] = 1;
	auto cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 4);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 49844044907.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 100000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.06, 0.01);
	const auto& reg = cluster.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::He] = 0;
	comp[Spec::V] = 2;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.17804, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 37523834.9, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 50000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.62, 0.01);
	const auto& regBis = cluster.getRegion();
	lo = regBis.getOrigin();
	hi = regBis.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::V] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.14131, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 1934083614.2, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 100000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.34, 0.01);
	const auto& regTer = cluster.getRegion();
	lo = regTer.getOrigin();
	hi = regTer.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::I] = 0;
	comp[Spec::He] = 2;
	comp[Spec::V] = 3;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 11);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.2038056, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& regQua = cluster.getRegion();
	lo = regQua.getOrigin();
	hi = regQua.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 3);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 4);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);
}

BOOST_AUTO_TEST_CASE(grouped)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 8 1" << std::endl
			  << "process=reaction" << std::endl
			  << "grouping=4 4 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char** argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	opts.readParams(argc, argv);

	using NetworkType = FeReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe = opts.getMaxImpurity();
	NetworkType::AmountType groupingWidthHe = opts.getGroupingWidthA();
	NetworkType::AmountType groupingWidthV = opts.getGroupingWidthB();
	// Adapt maxHe and maxV
	int i = 0;
	while (maxHe + 1 > pow(groupingWidthHe, i)) {
		++i;
	}
	maxHe = pow(groupingWidthV, i) - 1;
	i = 0;
	while (maxV + 1 > pow(groupingWidthV, i)) {
		++i;
	}
	maxV = pow(groupingWidthV, i) - 1;
	NetworkType network({maxHe, maxV, maxI},
		{{groupingWidthHe, groupingWidthV, maxI + 1}}, 1, opts);

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 27);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 33);
	// TODO: check it is within a given range?
	auto deviceMemorySize = network.getDeviceMemorySize();
	BOOST_REQUIRE(deviceMemorySize > 300000);
	BOOST_REQUIRE(deviceMemorySize < 330000);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 27);

	// Check the single vacancy
	auto vacancy = network.getSingleVacancy();
	BOOST_REQUIRE_EQUAL(vacancy.getId(), 0);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 3, 4, 27, 5, 6, 28, 7, 8, 29, 9, 10, 30, 11, 31, 12,
		32, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 2};
	knownDFill[1] = {1, 0, 2, 4, 27, 6, 28, 8, 29, 9, 10, 30, 11, 31, 12, 32,
		13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 3, 5, 7};
	knownDFill[2] = {2, 0, 26, 1, 3, 5, 7, 9, 13, 17, 21, 22, 4, 27};
	knownDFill[3] = {3, 0, 1, 18, 2, 4, 27, 5, 6, 28, 9, 10, 30, 11, 31, 13, 14,
		15, 17, 19, 21, 22, 23, 24, 26, 7, 8, 29, 12, 32, 16, 20, 25};
	knownDFill[4] = {4, 0, 27, 10, 30, 1, 2, 3, 5, 9, 13, 22, 6, 28};
	knownDFill[5] = {
		5, 0, 1, 19, 2, 3, 4, 27, 9, 10, 30, 13, 14, 17, 18, 21, 22, 23, 26, 7};
	knownDFill[6] = {6, 0, 28, 11, 31, 1, 2, 5, 3, 4, 27, 9, 13, 22, 8, 29};
	knownDFill[7] = {7, 0, 1, 20, 2, 3, 5, 9, 13, 17, 21, 22, 26};
	knownDFill[8] = {8, 0, 29, 12, 32, 1, 2, 7, 3, 6, 28, 4, 27, 5, 9, 13, 22};
	knownDFill[9] = {9, 0, 22, 1, 21, 3, 4, 27, 5, 6, 28, 7, 8, 29, 10, 30, 11,
		31, 12, 32, 13, 14, 15, 16, 18, 19, 20, 23, 24, 25};
	knownDFill[10] = {
		10, 0, 30, 23, 1, 4, 27, 3, 9, 17, 21, 26, 5, 14, 18, 13, 22, 11, 31};
	knownDFill[11] = {11, 0, 31, 24, 1, 6, 28, 3, 10, 30, 5, 9, 17, 21, 26, 15,
		19, 13, 22, 12, 32};
	knownDFill[12] = {12, 0, 32, 25, 1, 8, 29, 3, 11, 31, 5, 10, 30, 7, 9, 17,
		21, 26, 16, 20, 13, 22};
	knownDFill[13] = {13, 0, 1, 22, 3, 4, 27, 5, 6, 28, 7, 8, 29, 9, 10, 30, 11,
		31, 12, 32, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25};
	knownDFill[14] = {14, 0, 18, 1, 23, 3, 13, 5, 9, 22, 15};
	knownDFill[15] = {15, 0, 19, 1, 24, 3, 14, 5, 13, 9, 22, 16};
	knownDFill[16] = {16, 0, 20, 1, 25, 3, 15, 5, 14, 7, 13, 9, 22};
	knownDFill[17] = {17, 0, 21, 1, 26, 3, 5, 7, 9, 13, 22, 10, 30};
	knownDFill[18] = {18, 0, 3, 1, 14, 5, 9, 13, 22, 19};
	knownDFill[19] = {19, 0, 5, 1, 15, 3, 18, 9, 13, 22, 20};
	knownDFill[20] = {20, 0, 7, 1, 16, 3, 19, 5, 18, 9, 13, 22};
	knownDFill[21] = {21, 0, 9, 1, 17, 3, 5, 7, 13, 22, 10, 30};
	knownDFill[22] = {22, 0, 13, 1, 9, 3, 4, 27, 5, 6, 28, 7, 8, 29, 10, 30, 11,
		31, 12, 32, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25};
	knownDFill[23] = {23, 0, 14, 1, 10, 30, 3, 22, 5, 9, 13, 18, 24};
	knownDFill[24] = {24, 0, 15, 1, 11, 31, 3, 23, 5, 22, 9, 13, 19, 25};
	knownDFill[25] = {25, 0, 16, 1, 12, 32, 3, 24, 5, 23, 7, 22, 9, 13, 20};
	knownDFill[26] = {26, 0, 17, 1, 2, 3, 5, 7, 9, 22, 13, 21, 10, 30};
	knownDFill[27] = {27, 0, 4, 10, 30, 1, 2, 3, 5, 9, 13, 22, 6, 28};
	knownDFill[28] = {28, 0, 6, 11, 31, 1, 2, 5, 3, 4, 27, 9, 13, 22, 8, 29};
	knownDFill[29] = {29, 0, 8, 12, 32, 1, 2, 7, 3, 6, 28, 4, 27, 5, 9, 13, 22};
	knownDFill[30] = {
		30, 0, 10, 23, 1, 4, 27, 3, 9, 17, 21, 26, 5, 14, 18, 13, 22, 11, 31};
	knownDFill[31] = {31, 0, 11, 24, 1, 6, 28, 3, 10, 30, 5, 9, 17, 21, 26, 15,
		19, 13, 22, 12, 32};
	knownDFill[32] = {32, 0, 12, 25, 1, 8, 29, 3, 11, 31, 5, 10, 30, 7, 9, 17,
		21, 26, 16, 20, 13, 22};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 605);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	network.setTemperatures(temperatures);
	network.syncClusterDataOnHost();
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 751630489466.0, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 1), 48.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 3), 384.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::I, 1), 0.14131, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::He, 1), 90.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 540944269248, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-5.25014e+11, -5.24202e+11, -6.88497e+11,
		-1.54314e+13, -5.35907e+11, -3.97844e+12, 2.90214e+10, -3.42229e+11,
		6.08279e+11, -6.24652e+11, -1.85795e+11, 1.7113e+11, 6.35424e+11,
		-5.70605e+11, -1.58436e+11, 1.5447e+11, 5.60267e+11, -6.46997e+11,
		-1.47436e+11, 1.43382e+11, 5.18581e+11, -6.34698e+11, -6.48285e+11,
		-1.77674e+11, 1.51497e+11, 5.78836e+11, -6.64434e+11, -6.23897e+11,
		-2.23332e+10, 5.69677e+11, -5.35117e+11, -1.80248e+10, 5.0218e+11, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-5.39296e+11, 8.94081e+08,
		-2.76652e+11, -1.69531e+09, 8.47769e+08, -1.45181e+11, -1.69551e+09,
		8.47769e+08, -9.90247e+10, -1.69554e+09, 8.47769e+08, -4.40779e+08,
		-8.22511e+08, 0, -8.22635e+08, 0, -8.22637e+08, 0, 1.2078e+10,
		-1.68579e+08, -1.68588e+08, -1.68588e+08, -2.07355e+08, 7.84608e+08,
		-1.49194e+08, -1.49198e+08, -1.82928e+08, -1.99721e+09, -1.82172e+08,
		-1.8219e+08, -1.8219e+08, -2.15732e+08, 983384, -5.24202e+11,
		-7.01816e+09, -1.03034e+10, -8.92127e+10, 0, -8.92127e+10, 0,
		-8.92127e+10, 0, -9.32399e+09, -3.78735e+10, 0, -3.78735e+10, 0,
		-3.78735e+10, 0, -7.91224e+09, -7.76165e+09, -7.76165e+09, -7.76165e+09,
		-9.67534e+09, -6.86896e+09, -6.86896e+09, -6.86896e+09, -9.30735e+09,
		-1.03421e+10, -8.38785e+09, -8.38785e+09, -8.38785e+09, -1.00044e+10,
		1.10364e-10, 0.00140902, 38.3638, -6.92408e+11, 2.17302e+08,
		2.17302e+08, -1.03034e+10, -3.6493e+11, -1.89192e+11, -1.27981e+11,
		2.14787e+09, 2.05034e+08, 2.05034e+08, 2.52231e+09, 2.52231e+09,
		4.30157e+07, 0, -1.6943e+13, -2.76652e+11, 6.86896e+09, -2.68616e+11,
		-3.6493e+11, -3.09426e+12, 0, -5.83904e+11, -3.09431e+12, 0,
		-3.29038e+11, -1.37342e+12, 0, -1.37364e+12, 0, -2.9965e+11,
		-2.99407e+11, -2.99425e+11, -3.48743e+11, -2.76419e+11, -3.3926e+11,
		-3.18416e+11, -3.1553e+11, -3.15563e+11, -3.57224e+11, 4.50762e+08,
		1.65598e-11, 0, 7.35128e-12, 0, 1.60243e-12, 1.47931e-12, 1.6888e-12,
		-5.88309e+11, 6.42685e+07, 4.87882e+09, 1.31871e+09, 5.68783e+08,
		3.48487e+09, 4.56163e+10, -3.41172e+11, -2.00121e+11, 2.82814e+08,
		1.00509e+08, 1.4399e+09, 2.06997e-12, 0, -4.73739e+12, -1.45181e+11,
		6.86896e+09, 6.86896e+09, -1.89192e+11, 9.19357e+11, -1.60097e+12, 0,
		-1.71586e+11, -7.13727e+11, 0, -1.56675e+11, -1.56439e+11, -1.81099e+11,
		-1.44936e+11, -1.76357e+11, -1.67483e+11, -1.64509e+11, -1.85339e+11,
		4.50762e+08, -3.88182e+11, 6.42685e+07, 4.87882e+09, 1.31871e+09,
		5.68783e+08, 3.48487e+09, 2.3649e+10, 2.3649e+10, 0, 3.86788e+11, 0,
		2.82814e+08, 1.00509e+08, 1.4399e+09, 2.06997e-12, 0, -9.33002e+11,
		-9.90247e+10, 6.86896e+09, 6.86896e+09, -1.27981e+11, 5.83904e+11,
		5.83904e+11, -1.16601e+11, -1.06607e+11, -1.2264e+11, -1.1951e+11,
		-1.1475e+11, -1.25438e+11, -1.39395e+09, 6.42685e+07, 4.87882e+09,
		1.31871e+09, 5.68783e+08, 3.48487e+09, 1.59977e+10, 1.59977e+10,
		3.86788e+11, 3.86788e+11, 0, 2.00121e+11, 0, 2.00121e+11, 2.82814e+08,
		1.00509e+08, 1.4399e+09, -6.38573e+11, 1.50576e+09, -8.00248e+08,
		-1.66375e+07, 9.32659e+09, -3.29038e+11, -2.59408e+09, 5.18817e+09,
		-1.71586e+11, -2.59408e+09, 5.18817e+09, -1.16601e+11, -2.59408e+09,
		5.18817e+09, -2.21078e+09, 0, -2.26252e+09, 0, -2.26252e+09, 0, 316775,
		-4.8158e+08, -4.8158e+08, -4.8158e+08, -4.37618e+08, -4.37618e+08,
		-4.37618e+08, -5.12418e+08, -5.12418e+08, -5.12418e+08, -5.269e+11,
		-8.29896e+07, 2.41307e+09, 7.95687e+08, -5.7862e+09, 2.7879e+09,
		-9.75764e+09, 1.56888e+08, 8.20518e+10, 8.71858e+10, 8.4815e+10,
		8.93059e+10, -1.78432e+11, 7.02919e+08, 5.97957e+08, -1.14028e+08,
		-1.27368e+09, 1.83782e-12, 0, -3.48414e+11, -8.29896e+07, 2.41307e+09,
		7.95687e+08, -5.7862e+09, 2.7879e+09, -9.75764e+09, 0, 3.43409e+11, 0,
		1.78595e+11, 4.26887e+10, 4.52747e+10, 4.40892e+10, 4.63347e+10,
		7.02919e+08, 5.97957e+08, -1.14028e+08, -1.27368e+09, 1.83782e-12, 0,
		-5.00451e+09, -8.29896e+07, 2.41307e+09, 7.95687e+08, -5.7862e+09,
		2.7879e+09, -9.75764e+09, 3.43409e+11, 3.43409e+11, 0, 1.78432e+11,
		1.78432e+11, 0, 1.21047e+11, 2.89425e+10, 3.06599e+10, 2.98775e+10,
		3.13595e+10, 7.02919e+08, 5.97957e+08, -1.14028e+08, -1.27368e+09,
		-5.82074e+11, 2.77618e+08, 2.42983e+09, 8.139e+09, -2.9965e+11,
		-1.40205e+09, 1.40205e+09, -1.56675e+11, -1.40205e+09, 1.40205e+09,
		-1.06607e+11, -1.40205e+09, 1.40205e+09, -6.71307e+08, -8.04074e+08, 0,
		-8.04074e+08, 0, -8.04074e+08, 0, -1.50371e+08, -1.67906e+08,
		-1.67906e+08, -2.05034e+08, -1.50587e+08, -1.50587e+08, -1.50587e+08,
		-1.97895e+08, -1.80055e+08, -1.80055e+08, -1.80055e+08, -4.66624e+11,
		-1.93899e+07, 1.49198e+08, 6.26203e+08, 8.38787e+09, 2.25415e+08,
		2.99482e+11, -1.56439e+11, -4.8158e+08, -2.16219e+09, 1.60243e-12,
		-3.10167e+11, -1.93899e+07, 1.49198e+08, 6.26203e+08, 8.38785e+09, 0,
		2.99425e+11, 1.56675e+11, 1.56507e+11, -4.8158e+08, -2.16219e+09,
		1.60243e-12, -1.07419e+10, -1.93899e+07, 1.49198e+08, 6.26203e+08,
		8.38785e+09, 2.99425e+11, 2.99425e+11, 1.56439e+11, 1.56439e+11,
		1.06607e+11, 1.06439e+11, -4.8158e+08, -2.16219e+09, -6.62575e+11,
		-7.9929e+06, 2.02162e+08, 3.29059e+08, 1.0006e+10, -3.48743e+11,
		-1.81099e+11, -1.2264e+11, 6.71307e+08, 4.66273e+08, 9.23234e+09,
		5.48327e+07, 0, -4.3185e+11, 2.76503e+11, 2.32972e+08, 8.92693e+08,
		7.76166e+09, -1.44936e+11, -4.37618e+08, -1.50587e+08, -1.95421e+09,
		1.47931e-12, -2.8598e+11, 1.45032e+11, 1.45181e+11, 8.92693e+08,
		7.76165e+09, 0, 2.76419e+11, -4.37618e+08, -1.50587e+08, -1.95421e+09,
		1.47931e-12, -9.56057e+09, 9.88755e+10, 9.90247e+10, 8.92693e+08,
		7.76165e+09, 2.76419e+11, 2.76419e+11, 1.44936e+11, 1.44936e+11,
		-4.37618e+08, -1.50587e+08, -1.95421e+09, -6.47376e+11, 4.28474e+08,
		6.30636e+08, 3.67986e+08, 9.67814e+09, -3.3926e+11, -1.76357e+11,
		-1.1951e+11, 2.14435e+09, -1.80065e+08, 5.33416e+07, 0, -6.67384e+11,
		-1.81723e+09, -2.02307e+09, -1.01807e+09, 6.5772e+09, -3.18416e+11,
		-1.47599e+10, 2.21399e+10, -1.67483e+11, -1.47599e+10, 2.21399e+10,
		-1.1475e+11, -1.47599e+10, 2.21399e+10, -1.02393e+10, 0, -1.02393e+10,
		0, -1.02393e+10, 0, -2.16219e+09, -2.16219e+09, -2.16219e+09,
		-1.95421e+09, -1.95421e+09, -1.95421e+09, -2.52231e+09, -2.27513e+09,
		-2.30808e+09, -2.30808e+09, -4.91675e+11, -1.36016e+07, 1.68588e+08,
		-1.3122e+10, 9.4684e+09, -1.42026e+10, 2.85278e+09, 3.16108e+11,
		-1.64509e+11, -5.12418e+08, -2.94686e+07, 1.50587e+08, 1.6888e-12,
		-3.27134e+11, -1.36016e+07, 1.68588e+08, -1.3122e+10, 9.46837e+09,
		-1.42026e+10, 0, 3.15563e+11, 1.67483e+11, 1.65175e+11, -5.12418e+08,
		-2.94686e+07, 1.50587e+08, 1.6888e-12, -1.15706e+10, -1.36016e+07,
		1.68588e+08, -1.3122e+10, 9.46837e+09, -1.42026e+10, 3.15563e+11,
		3.15563e+11, 1.64509e+11, 1.64509e+11, 1.1475e+11, 1.12442e+11,
		-5.12418e+08, -2.94686e+07, 1.50587e+08, -6.78224e+11, -7.14739e+06,
		2.10155e+08, 2.99036e+08, 1.03044e+10, -3.57224e+11, -1.85339e+11,
		-1.25438e+11, 2.93665e+09, 2.93665e+09, 1.97895e+08, 1.97895e+08,
		5.6166e+07, 0, -5.93266e+11, -2.26607e+07, 4.82752e+08, -6.10705e+08,
		-9.13538e+07, -5.11114e+09, -3.04109e+10, -4.17199e+11, -2.00121e+11,
		-2.87458e+08, -5.50408e+07, -1.09091e+09, 0, 2.06997e-12, -3.93133e+11,
		-2.26607e+07, 4.7917e+08, -6.10705e+08, -9.13538e+07, -5.11114e+09,
		-1.5766e+10, -1.5766e+10, 0, 0, 3.86788e+11, -2.87458e+08, -5.50408e+07,
		-1.09091e+09, 0, 2.06997e-12, -6.34432e+09, -2.26607e+07, 4.7917e+08,
		-6.10705e+08, -9.13538e+07, -5.11114e+09, -1.06651e+10, -1.06651e+10,
		3.86788e+11, 0, 3.86788e+11, 0, 2.00121e+11, 2.00121e+11, -2.87458e+08,
		-5.50408e+07, -1.09091e+09, -5.33689e+11, -1.47204e+08, -1.93195e+09,
		3.11871e+08, -1.83055e+10, 3.34548e+09, -1.17092e+10, -3.34005e+11,
		-9.9119e+10, 3.48743e+10, -3.3926e+10, 1.07167e+11, -1.78432e+11,
		-2.18433e+08, -6.30025e+08, -2.49294e+08, -3.51547e+09, 0, 1.83782e-12,
		-3.55204e+11, -1.47204e+08, -1.93046e+09, 3.11871e+08, -1.83055e+10,
		3.34548e+09, -1.17092e+10, 0, 0, 3.43409e+11, 4.60016e+09, -5.18832e+10,
		1.81099e+10, -1.76357e+10, 5.56016e+10, -2.18433e+08, -6.30025e+08,
		-2.49294e+08, -3.51547e+09, 0, 1.83782e-12, -1.17942e+10, -1.47204e+08,
		-1.93046e+09, 3.11871e+08, -1.83055e+10, 3.34548e+09, -1.17092e+10,
		3.43409e+11, 0, 3.43409e+11, 1.78432e+11, 0, 1.78432e+11, 2.96406e+09,
		-3.53878e+10, 1.2264e+10, -1.1951e+10, 3.76314e+10, -2.18433e+08,
		-6.30025e+08, -2.49294e+08, -3.51547e+09};
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					knownPartials[startingIdx + j], 0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::He] = 1;
	auto cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 3);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 49844044907.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 100000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.06, 0.01);
	const auto& reg = cluster.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::He] = 0;
	comp[Spec::V] = 2;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 13);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.17804, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 37523834.9, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 50000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.62, 0.01);
	const auto& regBis = cluster.getRegion();
	lo = regBis.getOrigin();
	hi = regBis.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::V] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 1);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.14131, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 1934083614.2, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 100000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.34, 0.01);
	const auto& regTer = cluster.getRegion();
	lo = regTer.getOrigin();
	hi = regTer.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::I] = 0;
	comp[Spec::He] = 2;
	comp[Spec::V] = 6;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 11);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.248263, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& regQua = cluster.getRegion();
	lo = regQua.getOrigin();
	hi = regQua.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 4);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 8);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);
}

BOOST_AUTO_TEST_SUITE_END()
