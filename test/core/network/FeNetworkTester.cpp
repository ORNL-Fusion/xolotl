#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>

using namespace std;
using namespace xolotl;
using namespace core;
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
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=3 0 0 3 1" << std::endl
			  << "process=reaction sink" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

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
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 8 1" << std::endl
			  << "process=reaction" << std::endl
			  << "grouping=4 4 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

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

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 30);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 36);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 30);

	// Check the single vacancy
	auto vacancy = network.getSingleVacancy();
	BOOST_REQUIRE_EQUAL(vacancy.getId(), 0);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 3, 5, 7, 9, 10, 30, 11, 31, 12, 32, 13, 14, 15, 16,
		17, 18, 33, 19, 34, 20, 35, 21, 22, 23, 24, 25, 26, 27, 28, 29, 4, 6, 8,
		2};
	knownDFill[1] = {1, 0, 2, 4, 6, 8, 9, 10, 30, 11, 31, 12, 32, 13, 14, 15,
		16, 17, 18, 33, 19, 34, 20, 35, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3,
		5, 7};
	knownDFill[2] = {2, 0, 29, 1, 3, 5, 7, 9, 13, 17, 24, 25, 4};
	knownDFill[3] = {3, 0, 1, 21, 2, 4, 5, 6, 9, 10, 30, 11, 31, 13, 14, 15, 17,
		18, 33, 19, 34, 22, 24, 25, 26, 27, 29, 7, 8, 12, 32, 16, 20, 35, 23,
		28};
	knownDFill[4] = {4, 0, 18, 33, 1, 2, 3, 5, 9, 10, 30, 25, 13, 6};
	knownDFill[5] = {5, 0, 1, 22, 2, 3, 4, 9, 10, 30, 13, 14, 17, 18, 33, 21,
		24, 25, 26, 29, 7};
	knownDFill[6] = {6, 0, 19, 34, 1, 2, 5, 3, 4, 9, 11, 31, 25, 13, 8};
	knownDFill[7] = {7, 0, 1, 23, 2, 3, 5, 9, 13, 17, 24, 25, 29};
	knownDFill[8] = {8, 0, 20, 35, 1, 2, 7, 3, 6, 4, 5, 9, 12, 32, 25, 13};
	knownDFill[9] = {9, 0, 25, 1, 24, 3, 5, 7, 10, 30, 11, 31, 12, 32, 13, 14,
		15, 16, 21, 22, 23, 26, 27, 28};
	knownDFill[10] = {
		10, 0, 30, 26, 1, 18, 33, 3, 9, 24, 5, 21, 13, 25, 14, 11, 31};
	knownDFill[11] = {
		11, 0, 31, 27, 1, 19, 34, 3, 10, 30, 5, 9, 24, 22, 13, 25, 15, 12, 32};
	knownDFill[12] = {12, 0, 32, 28, 1, 20, 35, 3, 11, 31, 5, 10, 30, 7, 9, 24,
		23, 13, 25, 16};
	knownDFill[13] = {13, 0, 1, 25, 3, 5, 7, 9, 10, 30, 11, 31, 12, 32, 14, 15,
		16, 17, 18, 33, 19, 34, 20, 35, 21, 22, 23, 24, 26, 27, 28};
	knownDFill[14] = {14, 0, 21, 1, 26, 3, 13, 5, 9, 25, 15};
	knownDFill[15] = {15, 0, 22, 1, 27, 3, 14, 5, 13, 9, 25, 16};
	knownDFill[16] = {16, 0, 23, 1, 28, 3, 15, 5, 14, 7, 13, 9, 25};
	knownDFill[17] = {17, 0, 24, 1, 29, 3, 5, 7, 9, 13, 25, 18, 33};
	knownDFill[18] = {
		18, 0, 10, 30, 33, 1, 4, 3, 17, 29, 5, 9, 14, 26, 13, 25, 19, 34};
	knownDFill[19] = {19, 0, 11, 31, 34, 1, 6, 3, 18, 33, 5, 17, 29, 9, 15, 27,
		13, 25, 20, 35};
	knownDFill[20] = {20, 0, 12, 32, 35, 1, 8, 3, 19, 34, 5, 18, 33, 7, 17, 29,
		9, 16, 28, 13, 25};
	knownDFill[21] = {21, 0, 3, 1, 14, 5, 9, 13, 25, 22};
	knownDFill[22] = {22, 0, 5, 1, 15, 3, 21, 9, 13, 25, 23};
	knownDFill[23] = {23, 0, 7, 1, 16, 3, 22, 5, 21, 9, 13, 25};
	knownDFill[24] = {24, 0, 9, 1, 17, 3, 5, 7, 13, 25, 10, 30};
	knownDFill[25] = {25, 0, 13, 1, 9, 3, 5, 7, 10, 30, 11, 31, 12, 32, 14, 15,
		16, 21, 22, 23, 24, 26, 27, 28};
	knownDFill[26] = {26, 0, 14, 1, 10, 30, 3, 25, 5, 9, 13, 21, 27};
	knownDFill[27] = {27, 0, 15, 1, 11, 31, 3, 26, 5, 25, 9, 13, 22, 28};
	knownDFill[28] = {28, 0, 16, 1, 12, 32, 3, 27, 5, 26, 7, 25, 9, 13, 23};
	knownDFill[29] = {29, 0, 17, 1, 2, 3, 5, 7, 9, 25, 13, 24, 18, 33};
	knownDFill[30] = {
		30, 0, 10, 26, 1, 18, 33, 3, 9, 24, 5, 21, 13, 25, 14, 11, 31};
	knownDFill[31] = {
		31, 0, 11, 27, 1, 19, 34, 3, 10, 30, 5, 9, 24, 22, 13, 25, 15, 12, 32};
	knownDFill[32] = {32, 0, 12, 28, 1, 20, 35, 3, 11, 31, 5, 10, 30, 7, 9, 24,
		23, 13, 25, 16};
	knownDFill[33] = {
		33, 0, 10, 30, 18, 1, 4, 3, 17, 29, 5, 9, 14, 26, 13, 25, 19, 34};
	knownDFill[34] = {34, 0, 11, 31, 19, 1, 6, 3, 18, 33, 5, 17, 29, 9, 15, 27,
		13, 25, 20, 35};
	knownDFill[35] = {35, 0, 12, 32, 20, 1, 8, 3, 19, 34, 5, 18, 33, 7, 17, 29,
		9, 16, 28, 13, 25};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 659);
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
		network.getTotalConcentration(dConcs, Spec::He, 1), 27.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 3), 132.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::I, 1), 0.14131, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::He, 1), 48.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 535857652671, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-5.22471e+11, -2.87474e+11, -6.88335e+11,
		-9.97246e+12, -1.95141e+11, -2.56667e+12, 1.83297e+11, -3.42229e+11,
		6.76208e+11, -6.26418e+11, -1.73605e+11, 1.73747e+11, 6.25311e+11,
		-5.69668e+11, -1.58436e+11, 1.5447e+11, 5.60267e+11, -6.47052e+11,
		-1.78184e+11, 1.8831e+11, 6.65332e+11, -1.47436e+11, 1.43382e+11,
		5.18581e+11, -6.3465e+11, -6.54547e+11, -1.68392e+11, 1.6078e+11,
		5.88118e+11, -6.6449e+11, -5.17148e+11, -1.50675e+10, 4.9055e+11,
		-5.43221e+11, -1.1202e+10, 5.23559e+11, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-5.36753e+11, 8.94081e+08,
		-2.76652e+11, -1.45181e+11, -9.90247e+10, -4.40779e+08, -3.95124e+08,
		1135.1, -3.95179e+08, 8.91754, -3.9518e+08, 0.140363, 1.2078e+10,
		-1.68579e+08, -1.68588e+08, -1.68588e+08, -2.07355e+08, -4.27392e+08,
		1269.95, -4.27455e+08, 33.6357, -4.27457e+08, 1.33481, 7.84608e+08,
		-1.49194e+08, -1.49198e+08, -1.82928e+08, -1.99721e+09, -1.82172e+08,
		-1.8219e+08, -1.8219e+08, -2.15732e+08, 31733.7, 1492.92, 94.0114,
		983384, -2.87473e+11, -7.01816e+09, -1.03034e+10, -1.03034e+10,
		-1.03034e+10, -1.03034e+10, -9.32399e+09, -1.81937e+10, 1.98817e-15,
		-1.81937e+10, 1.98817e-15, -1.81937e+10, 7.09636e-14, -7.91224e+09,
		-7.76165e+09, -7.76165e+09, -7.76165e+09, -9.67534e+09, -1.96797e+10,
		1.24043e-15, -1.96797e+10, 1.24043e-15, -1.96797e+10, 4.87899e-15,
		-6.86896e+09, -6.86896e+09, -6.86896e+09, -9.30735e+09, -1.03421e+10,
		-8.38785e+09, -8.38785e+09, -8.38785e+09, -1.00044e+10, 1.10364e-10,
		0.00140902, 38.3638, -6.92408e+11, 2.17302e+08, 2.17302e+08,
		-1.03034e+10, -3.6493e+11, -1.89192e+11, -1.27981e+11, 2.14787e+09,
		2.05034e+08, 2.05034e+08, 2.52231e+09, 2.52231e+09, 5.41836e+07,
		-1.14843e+13, -2.76652e+11, 6.86896e+09, -2.68616e+11, -3.6493e+11,
		-3.64876e+11, -5.83904e+11, -3.6493e+11, -3.29038e+11, -6.67572e+11,
		753412, -6.6767e+11, 0, -2.9965e+11, -2.99407e+11, -2.99425e+11,
		-3.48743e+11, -7.05856e+11, 667319, -7.05967e+11, 0, -2.76419e+11,
		-3.3926e+11, -3.18416e+11, -3.1553e+11, -3.15563e+11, -3.57224e+11,
		4.50762e+08, 1.95299e-12, 3.57316e-12, 0, 1.60243e-12, 3.77812e-12, 0,
		1.47931e-12, 1.6888e-12, -5.6448e+11, 3.20593e+08, 4.21955e+08,
		2.75113e+06, -1.03034e+10, 3.6493e+11, 0, -1.89192e+11, 2.73667e+08,
		3.0206e+09, 9.62969e+08, 3.70991e+09, 1.04113e+08, 1.95299e-12,
		-3.32562e+12, -1.45181e+11, 6.86896e+09, 6.86896e+09, -1.89192e+11,
		9.19357e+11, -1.89192e+11, -1.71586e+11, -3.47289e+11, 0, -1.56675e+11,
		-1.56439e+11, -1.81099e+11, -3.66437e+11, 0, -1.44936e+11, -1.76357e+11,
		-1.67483e+11, -1.64509e+11, -1.85339e+11, 4.50762e+08, -3.75234e+11,
		3.20593e+08, 4.21955e+08, 2.75113e+06, -1.03034e+10, 1.89192e+11,
		1.89192e+11, 0, 3.6493e+11, 2.73667e+08, 3.0206e+09, 9.62969e+08,
		3.70991e+09, 1.04113e+08, 1.95299e-12, -9.33002e+11, -9.90247e+10,
		6.86896e+09, 6.86896e+09, -1.27981e+11, 5.83904e+11, 5.83904e+11,
		-1.16601e+11, -1.06607e+11, -1.2264e+11, -1.1951e+11, -1.1475e+11,
		-1.25438e+11, -1.03034e+10, 3.20593e+08, 4.21955e+08, 2.75113e+06,
		-1.03034e+10, 1.27981e+11, 1.27981e+11, 3.6493e+11, 3.6493e+11,
		1.89192e+11, 1.89192e+11, 2.73667e+08, 3.0206e+09, 9.62969e+08,
		3.70991e+09, 1.04113e+08, -6.40389e+11, 1.50576e+09, -8.00248e+08,
		-1.66375e+07, 9.32659e+09, -3.29038e+11, -1.71586e+11, -1.16601e+11,
		-4.98829e+08, 2.49415e+08, -5.47334e+08, 2.73667e+08, -5.47334e+08,
		2.73667e+08, 316775, -4.8158e+08, -4.8158e+08, -4.8158e+08,
		-4.37618e+08, -4.37618e+08, -4.37618e+08, -5.12418e+08, -5.12418e+08,
		-5.12418e+08, -5.15117e+11, -5.70976e+07, 2.36128e+09, 1.81122e+08,
		1.85752e+08, 4.91995e+09, -2.45998e+09, 3.13776e+08, 1.64601e+11,
		1.6963e+11, -1.73645e+11, 1.19591e+09, -1.98303e+07, -4.1507e+08,
		1.16505e+09, 1.78658e-12, 0, -3.41423e+11, -5.70976e+07, 2.36165e+09,
		1.81122e+08, 1.85752e+08, 4.91994e+09, -2.45997e+09, 0, 3.33835e+11, 0,
		1.73971e+11, 8.58747e+10, 8.81785e+10, 1.19591e+09, -1.98303e+07,
		-4.1507e+08, 1.16505e+09, 1.78658e-12, 0, -7.58798e+09, -5.70976e+07,
		2.36165e+09, 1.81122e+08, 1.85752e+08, 4.91993e+09, -2.45997e+09,
		3.33835e+11, 3.33835e+11, 0, 1.73645e+11, 1.73645e+11, 0, 1.18056e+11,
		5.83825e+10, 5.9755e+10, 1.19591e+09, -1.98303e+07, -4.1507e+08,
		1.16505e+09, -5.81137e+11, 2.77618e+08, 2.42983e+09, 8.139e+09,
		-2.9965e+11, -1.56675e+11, -1.06607e+11, -6.71307e+08, -3.87622e+08, 0,
		-3.87622e+08, 0, -3.87622e+08, 0, -1.50371e+08, -1.67906e+08,
		-1.67906e+08, -2.05034e+08, -2.08226e+08, 1.04113e+08, -2.08226e+08,
		1.04113e+08, -2.08226e+08, 1.04113e+08, -1.50587e+08, -1.50587e+08,
		-1.50587e+08, -1.97895e+08, -1.80055e+08, -1.80055e+08, -1.80055e+08,
		-4.66624e+11, -1.93899e+07, 1.49198e+08, 6.26203e+08, 8.38787e+09,
		2.25415e+08, 2.99482e+11, -1.56439e+11, -4.8158e+08, -2.16219e+09,
		1.60243e-12, -3.10167e+11, -1.93899e+07, 1.49198e+08, 6.26203e+08,
		8.38785e+09, 0, 2.99425e+11, 1.56675e+11, 1.56507e+11, -4.8158e+08,
		-2.16219e+09, 1.60243e-12, -1.07419e+10, -1.93899e+07, 1.49198e+08,
		6.26203e+08, 8.38785e+09, 2.99425e+11, 2.99425e+11, 1.56439e+11,
		1.56439e+11, 1.06607e+11, 1.06439e+11, -4.8158e+08, -2.16219e+09,
		-6.62575e+11, -7.9929e+06, 2.02162e+08, 3.29059e+08, 1.0006e+10,
		-3.48743e+11, -1.81099e+11, -1.2264e+11, 6.71307e+08, 4.66273e+08,
		9.23234e+09, 5.48854e+07, -2.74427e+07, -5.41389e+11, -1.21039e+07,
		1.52924e+09, -5.6892e+08, 2.45827e+09, 2.69175e+09, 5.15173e+09,
		-3.05176e-05, 1.74372e+11, 1.78612e+11, -1.83219e+11, 4.96999e+08,
		2.4079e+08, 1.41025e+09, 1.41754e+08, 1.77236e+09, 1.88906e-12, 0,
		-3.58114e+11, -1.21039e+07, 1.52924e+09, -5.6892e+08, 2.45859e+09,
		2.69175e+09, 5.15172e+09, 0, 3.52984e+11, 0, 1.83219e+11, 9.05493e+10,
		9.26694e+10, 4.96999e+08, 2.4079e+08, 1.41025e+09, 1.41754e+08,
		1.77236e+09, 1.88906e-12, 0, -5.13091e+09, -1.21039e+07, 1.52924e+09,
		-5.6892e+08, 2.45859e+09, 2.69175e+09, 5.15172e+09, 3.52984e+11,
		3.52984e+11, 0, 1.83219e+11, 1.83219e+11, 0, 1.24039e+11, 6.13198e+10,
		6.2719e+10, 4.96999e+08, 2.4079e+08, 1.41025e+09, 1.41754e+08,
		1.77236e+09, -4.3185e+11, 2.76503e+11, 2.32972e+08, 8.92693e+08,
		7.76166e+09, -1.44936e+11, -4.37618e+08, -1.50587e+08, -1.95421e+09,
		1.47931e-12, -2.8598e+11, 1.45032e+11, 1.45181e+11, 8.92693e+08,
		7.76165e+09, 0, 2.76419e+11, -4.37618e+08, -1.50587e+08, -1.95421e+09,
		1.47931e-12, -9.56057e+09, 9.88755e+10, 9.90247e+10, 8.92693e+08,
		7.76165e+09, 2.76419e+11, 2.76419e+11, 1.44936e+11, 1.44936e+11,
		-4.37618e+08, -1.50587e+08, -1.95421e+09, -6.47376e+11, 4.28474e+08,
		6.30636e+08, 3.67986e+08, 9.67814e+09, -3.3926e+11, -1.76357e+11,
		-1.1951e+11, 2.14435e+09, -1.80065e+08, 5.00117e+07, 2.50058e+07,
		-6.73645e+11, -1.81723e+09, -2.02307e+09, -1.01807e+09, 6.5772e+09,
		-3.18416e+11, -1.67483e+11, -1.1475e+11, -4.94654e+09, 0, -4.94654e+09,
		0, -4.94654e+09, 0, -2.16219e+09, -2.16219e+09, -2.16219e+09,
		-1.95421e+09, -1.95421e+09, -1.95421e+09, -2.52231e+09, -2.27513e+09,
		-2.30808e+09, -2.30808e+09, -4.91675e+11, -1.36016e+07, 1.68588e+08,
		-3.83942e+09, 9.09689e+09, -4.54845e+09, 2.85278e+09, 3.16108e+11,
		-1.64509e+11, -5.12418e+08, -2.94686e+07, 1.50587e+08, 1.6888e-12,
		-3.27134e+11, -1.36016e+07, 1.68588e+08, -3.83942e+09, 9.09686e+09,
		-4.54843e+09, 0, 3.15563e+11, 1.67483e+11, 1.65175e+11, -5.12418e+08,
		-2.94686e+07, 1.50587e+08, 1.6888e-12, -1.15706e+10, -1.36016e+07,
		1.68588e+08, -3.83942e+09, 9.09686e+09, -4.54843e+09, 3.15563e+11,
		3.15563e+11, 1.64509e+11, 1.64509e+11, 1.1475e+11, 1.12442e+11,
		-5.12418e+08, -2.94686e+07, 1.50587e+08, -6.78224e+11, -7.14739e+06,
		2.10155e+08, 2.99036e+08, 1.03044e+10, -3.57224e+11, -1.85339e+11,
		-1.25438e+11, 2.93665e+09, 2.93665e+09, 1.97895e+08, 1.97895e+08,
		5.622e+07, 2.811e+07, -5.24412e+11, -2.80985e+08, -8.35348e+09,
		-2.13458e+06, -1.78222e+10, 9.8399e+09, -4.91995e+09, -3.23614e+11,
		-3.28327e+11, 3.3926e+11, -1.73645e+11, -1.51659e+09, -1.81662e+08,
		-2.26529e+09, 1.99429e+09, 0, 1.78658e-12, -3.50718e+11, -2.80985e+08,
		-8.35194e+09, -2.13458e+06, -1.78222e+10, 9.83987e+09, -4.91994e+09, 0,
		0, 3.33835e+11, 4.77137e+09, -1.70874e+11, 1.76357e+11, -1.51659e+09,
		-1.81662e+08, -2.26529e+09, 1.99429e+09, 0, 1.78658e-12, -1.68824e+10,
		-2.80985e+08, -8.35194e+09, -2.13458e+06, -1.78222e+10, 9.83987e+09,
		-4.91993e+09, 3.33835e+11, 0, 3.33835e+11, 1.73645e+11, 0, 1.73645e+11,
		2.90896e+09, -1.1589e+11, 1.1951e+11, -1.51659e+09, -1.81662e+08,
		-2.26529e+09, 1.99429e+09, -5.51442e+11, -4.03249e+08, 2.27568e+09,
		-1.14162e+09, -9.41929e+09, -1.42962e+10, 1.03035e+10, -3.44503e+11,
		-3.48743e+11, 3.57224e+11, -1.83219e+11, 3.0838e+07, -4.8158e+08,
		-1.79567e+09, 2.97924e+08, -1.07145e+09, 0, 1.88906e-12, -3.68168e+11,
		-4.03249e+08, 2.27568e+09, -1.14162e+09, -9.41792e+09, -1.42962e+10,
		1.03034e+10, 0, 0, 3.52984e+11, 4.24016e+09, -1.81099e+11, 1.85339e+11,
		3.0838e+07, -4.8158e+08, -1.79567e+09, 2.97924e+08, -1.07145e+09, 0,
		1.88906e-12, -1.51845e+10, -4.03249e+08, 2.27568e+09, -1.14162e+09,
		-9.41792e+09, -1.42962e+10, 1.03034e+10, 3.52984e+11, 0, 3.52984e+11,
		1.83219e+11, 0, 1.83219e+11, 2.79851e+09, -1.2264e+11, 1.25438e+11,
		3.0838e+07, -4.8158e+08, -1.79567e+09, 2.97924e+08, -1.07145e+09};
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
				XOLOTL_REQUIRE_CLOSE_ZT(hPartials[startingIdx + index],
					knownPartials[startingIdx + j], 0.01, 1.0e-4);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 19);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.2635485, 0.01);
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
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 6);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 8);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);
}

BOOST_AUTO_TEST_SUITE_END()
