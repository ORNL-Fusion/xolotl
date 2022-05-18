#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the NE network.
 */
BOOST_AUTO_TEST_SUITE(NENetwork_testSuite)

BOOST_AUTO_TEST_CASE(fullyRefined)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=20 0 0 0 0" << std::endl
			  << "process=reaction" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;
	NetworkType network(
		{(NetworkType::AmountType)opts.getMaxImpurity()}, 1, opts);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 20);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 20);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.547, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0409168, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
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
	knownDFill[0] = {
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	knownDFill[1] = {1, 0, 2};
	knownDFill[2] = {2, 0, 1, 3};
	knownDFill[3] = {3, 0, 2, 4};
	knownDFill[4] = {4, 0, 3, 5};
	knownDFill[5] = {5, 0, 4, 6};
	knownDFill[6] = {6, 0, 5, 7};
	knownDFill[7] = {7, 0, 6, 8};
	knownDFill[8] = {8, 0, 7, 9};
	knownDFill[9] = {9, 0, 8, 10};
	knownDFill[10] = {10, 0, 9, 11};
	knownDFill[11] = {11, 0, 10, 12};
	knownDFill[12] = {12, 0, 11, 13};
	knownDFill[13] = {13, 0, 12, 14};
	knownDFill[14] = {14, 0, 13, 15};
	knownDFill[15] = {15, 0, 14, 16};
	knownDFill[16] = {16, 0, 15, 17};
	knownDFill[17] = {17, 0, 16, 18};
	knownDFill[18] = {18, 0, 17, 19};
	knownDFill[19] = {19, 0, 18};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 94);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
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

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::Xe, 1), 20.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::Xe, 3), 207.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1), 11.9849,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Xe, 1), 0.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 7.98584e-5, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-8.5332e-05, 2.45933e-06, -2.38186e-07,
		-1.8962e-07, -1.60127e-07, -1.39969e-07, -1.25163e-07, -1.13743e-07,
		-1.04618e-07, -9.71289e-08, -9.08512e-08, -8.54991e-08, -8.08717e-08,
		-7.68238e-08, -7.32473e-08, -7.006e-08, -6.71984e-08, -6.46122e-08,
		-6.22614e-08, 4.8543e-06, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-9.62794e-05, -3.01432e-06,
		-3.25251e-06, -3.44213e-06, -3.60226e-06, -3.74222e-06, -3.86739e-06,
		-3.98113e-06, -4.08575e-06, -4.18288e-06, -4.27373e-06, -4.35923e-06,
		-4.4401e-06, -4.51692e-06, -4.59017e-06, -4.66023e-06, -4.72743e-06,
		-4.79204e-06, -4.8543e-06, 4.04978e-32, -3.01432e-06, 7.93299e-06,
		6.13381e-15, -3.25251e-06, -2.38186e-07, 3.01432e-06, 3.6376e-16,
		-3.44213e-06, -1.8962e-07, 3.25251e-06, 6.75237e-17, -3.60226e-06,
		-1.60127e-07, 3.44213e-06, 3.95562e-17, -3.74222e-06, -1.39969e-07,
		3.60226e-06, 1.26426e-18, -3.86739e-06, -1.25163e-07, 3.74222e-06,
		4.09399e-19, -3.98113e-06, -1.13743e-07, 3.86739e-06, 3.98904e-22,
		-4.08575e-06, -1.04618e-07, 3.98113e-06, 2.37715e-20, -4.18288e-06,
		-9.71289e-08, 4.08575e-06, 4.26867e-21, -4.27373e-06, -9.08512e-08,
		4.18288e-06, 1.31745e-23, -4.35923e-06, -8.54991e-08, 4.27373e-06,
		7.38572e-25, -4.4401e-06, -8.08717e-08, 4.35923e-06, 2.35721e-25,
		-4.51692e-06, -7.68238e-08, 4.4401e-06, 4.20612e-26, -4.59017e-06,
		-7.32473e-08, 4.51692e-06, 2.34921e-27, -4.66023e-06, -7.006e-08,
		4.59017e-06, 7.47349e-28, -4.72743e-06, -6.71984e-08, 4.66023e-06,
		7.30852e-30, -4.79204e-06, -6.46122e-08, 4.72743e-06, 4.07174e-31,
		-4.8543e-06, -6.22614e-08, 4.79204e-06, 4.04978e-32, -4.04978e-32,
		4.8543e-06, 4.8543e-06};
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
	comp[Spec::Xe] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 3.6298e-7, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.0, 0.01);
	const auto& reg = cluster.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 2);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Xe] = 10;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 9);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.61702, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 46.9, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& regBis = cluster.getRegion();
	lo = regBis.getOrigin();
	hi = regBis.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 10);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 11);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
}

BOOST_AUTO_TEST_CASE(grouped)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=25 0 0 0 0" << std::endl
			  << "grouping=11 4" << std::endl
			  << "process=reaction" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the boundaries from the options
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
	// Adapt maxXe
	int i = 0;
	while (maxXe + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxXe = pow(groupingWidth, i) - 1;
	NetworkType network({maxXe}, {{groupingWidth}}, 1, opts);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 16);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 19);

	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 16);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {
		0, 1, 16, 2, 3, 4, 17, 5, 18, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	knownDFill[1] = {1, 0, 16, 4, 17, 5, 18};
	knownDFill[2] = {2, 0, 8, 9};
	knownDFill[3] = {3, 0, 11, 12};
	knownDFill[4] = {4, 0, 17, 14, 1, 16};
	knownDFill[5] = {5, 0, 1, 16, 18, 6};
	knownDFill[6] = {6, 0, 5, 18, 15};
	knownDFill[7] = {7, 0, 8};
	knownDFill[8] = {8, 0, 7, 2};
	knownDFill[9] = {9, 0, 2, 10};
	knownDFill[10] = {10, 0, 9, 11};
	knownDFill[11] = {11, 0, 10, 3};
	knownDFill[12] = {12, 0, 3, 13};
	knownDFill[13] = {13, 0, 12, 14};
	knownDFill[14] = {14, 0, 13, 4, 17};
	knownDFill[15] = {15, 0, 6};
	knownDFill[16] = {16, 0, 1, 4, 17, 5, 18};
	knownDFill[17] = {17, 0, 4, 14, 1, 16};
	knownDFill[18] = {18, 0, 1, 16, 5, 6};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 105);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 4.15126e+103, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::Xe, 1), 25.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::Xe, 3), 322.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1), 16.0555,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Xe, 1), 0.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 0.00010499, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-2.07563e+103, -1.76185e-07, -1.8962e-07,
		-1.13743e-07, 5.18908e+102, -1.5084e-07, 7.36301e-06, 2.45933e-06,
		-2.38186e-07, -1.60127e-07, -1.39969e-07, -1.25163e-07, -1.04618e-07,
		-9.71289e-08, -2.07563e+103, 5.1366e-06, -5.49878e-06, -6.22689e+102,
		-5.81881e-06, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-0.000121415, -1.9034e-05, 5.44285e-41,
		-3.44213e-06, -3.98113e-06, 4.15126e+103, -6.22689e+103, -1.99994e-05,
		4.65987e-41, -5.1366e-06, -3.01432e-06, -3.25251e-06, -3.60226e-06,
		-3.74222e-06, -3.86739e-06, -4.08575e-06, -4.18288e-06, -4.27373e-06,
		2.2293e-36, -1.18963e-06, -1.76185e-07, -1.78444e-06, 1.11915e-06,
		1.67873e-06, 1.53129e-40, -2.29693e-40, -3.44213e-06, -1.8962e-07,
		3.25251e-06, 6.75237e-17, -3.98113e-06, -1.13743e-07, 3.86739e-06,
		3.98904e-22, -1.03782e+103, -1.72945e-06, 1.55672e+103, 1.06843e-06,
		1.44058e-40, -2.16086e-40, -1.24996e-06, -1.5084e-07, 1.18963e-06,
		1.78444e-06, -1.87494e-06, 2.43742e-49, -5.1366e-06, 7.36301e-06,
		4.99984e-06, 7.49977e-06, 2.2293e-36, -3.01432e-06, 7.93299e-06,
		6.13381e-15, -3.25251e-06, -2.38186e-07, 3.01432e-06, 3.6376e-16,
		-3.60226e-06, -1.60127e-07, 3.44213e-06, 3.95562e-17, -3.74222e-06,
		-1.39969e-07, 3.60226e-06, 1.26426e-18, -3.86739e-06, -1.25163e-07,
		3.74222e-06, 4.09399e-19, -4.08575e-06, -1.04618e-07, 3.98113e-06,
		2.37715e-20, -4.18288e-06, -9.71289e-08, 4.08575e-06, 4.26867e-21,
		-4.27373e-06, -9.08512e-08, 4.18288e-06, 4.15126e+103, -6.22689e+103,
		-2.2293e-36, 5.1366e-06, 5.1366e-06, -3.56888e-06, -5.49878e-06,
		1.42755e-06, -1.34298e-06, -2.01447e-06, 1.83755e-40, -2.75632e-40,
		-1.86807e+103, -3.29659e-06, 1.24538e+103, -1.28212e-06, 1.72869e-40,
		-2.59304e-40, -3.74988e-06, -5.81881e-06, -1.42755e-06, -2.14133e-06,
		1.49995e-06, 2.92491e-49};
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
				XOLOTL_REQUIRE_CLOSE(hPartials[startingIdx + index],
					knownPartials[startingIdx + j], 0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::Xe] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 3.6298e-7, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.0, 0.01);
	const auto& reg = cluster.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 2);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Xe] = 20;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 5);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.7961, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 79.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& regBis = cluster.getRegion();
	lo = regBis.getOrigin();
	hi = regBis.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 20);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 24);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
}

BOOST_AUTO_TEST_CASE(fullyRefined_ReSo)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=20 0 0 0 0" << std::endl
			  << "process=reaction resolution" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;
	NetworkType network(
		{(NetworkType::AmountType)opts.getMaxImpurity()}, 1, opts);
	network.setFissionRate(8.0e-9);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == true);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	knownDFill[1] = {1, 0, 2};
	knownDFill[2] = {2, 0, 1, 3};
	knownDFill[3] = {3, 0, 2, 4};
	knownDFill[4] = {4, 0, 3, 5};
	knownDFill[5] = {5, 0, 4, 6};
	knownDFill[6] = {6, 0, 5, 7};
	knownDFill[7] = {7, 0, 6, 8};
	knownDFill[8] = {8, 0, 7, 9};
	knownDFill[9] = {9, 0, 8, 10};
	knownDFill[10] = {10, 0, 9, 11};
	knownDFill[11] = {11, 0, 10, 12};
	knownDFill[12] = {12, 0, 11, 13};
	knownDFill[13] = {13, 0, 12, 14};
	knownDFill[14] = {14, 0, 13, 15};
	knownDFill[15] = {15, 0, 14, 16};
	knownDFill[16] = {16, 0, 15, 17};
	knownDFill[17] = {17, 0, 16, 18};
	knownDFill[18] = {18, 0, 17, 19};
	knownDFill[19] = {19, 0, 18};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 94);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 0.117914, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 1.720317, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-1.83318, 0.0529346, -0.00516169,
		-0.00410766, -0.00346739, -0.00302979, -0.00270841, -0.00246059,
		-0.00226261, -0.00210017, -0.00196403, -0.001848, -0.0017477,
		-0.00165998, -0.00158248, -0.00151344, -0.00145146, -0.00139546,
		-0.00134456, 0.104389, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-2.07406, -0.064144, -0.069715,
		-0.0738305, -0.0773028, -0.0803359, -0.0830468, -0.0855092, -0.0877732,
		-0.0898745, -0.0918394, -0.0936881, -0.0954365, -0.0970969, -0.0986799,
		-0.100194, -0.101645, -0.103041, -0.104386, 0.000182462, -0.0653302,
		0.170893, 0.000350842, -0.0704167, -0.00513103, 0.0649348, 0.000320188,
		-0.0744708, -0.00408481, 0.0700659, 0.000297336, -0.0778975,
		-0.00344947, 0.0741507, 0.000279413, -0.0808948, -0.00301523, 0.0776001,
		0.000264854, -0.0835765, -0.00269627, 0.0806154, 0.000252718,
		-0.0860146, -0.00245026, 0.0833116, 0.000242396, -0.088258, -0.00225369,
		0.0857619, 0.000233478, -0.0903414, -0.00209236, 0.0880156, 0.000225672,
		-0.0922907, -0.00195713, 0.0901079, 0.000218764, -0.0941257,
		-0.00184183, 0.0920651, 0.000212596, -0.0958616, -0.00174215, 0.0939069,
		0.000207045, -0.097511, -0.00165494, 0.095649, 0.000202015, -0.0990839,
		-0.0015779, 0.097304, 0.000197429, -0.100589, -0.00150924, 0.0988819,
		0.000193228, -0.102032, -0.00144759, 0.100391, 0.000189359, -0.10342,
		-0.00139188, 0.101839, 0.000185782, -0.104758, -0.00134124, 0.103231,
		0.000182462, -0.000182462, 0.104572, 0.104572};
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
}

BOOST_AUTO_TEST_CASE(grouped_ReSo)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=25 0 0 0 0" << std::endl
			  << "grouping=11 4" << std::endl
			  << "process=reaction resolution" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the boundaries from the options
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
	// Adapt maxXe
	int i = 0;
	while (maxXe + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxXe = pow(groupingWidth, i) - 1;
	NetworkType network({maxXe}, {{groupingWidth}}, 1, opts);
	network.setFissionRate(8.0e-9);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {
		0, 1, 16, 2, 3, 4, 17, 5, 18, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	knownDFill[1] = {1, 0, 16, 4, 17, 5, 18};
	knownDFill[2] = {2, 0, 8, 9};
	knownDFill[3] = {3, 0, 11, 12};
	knownDFill[4] = {4, 0, 17, 14, 1, 16};
	knownDFill[5] = {5, 0, 1, 16, 18, 6};
	knownDFill[6] = {6, 0, 5, 18, 15};
	knownDFill[7] = {7, 0, 8};
	knownDFill[8] = {8, 0, 7, 2};
	knownDFill[9] = {9, 0, 2, 10};
	knownDFill[10] = {10, 0, 9, 11};
	knownDFill[11] = {11, 0, 10, 3};
	knownDFill[12] = {12, 0, 3, 13};
	knownDFill[13] = {13, 0, 12, 14};
	knownDFill[14] = {14, 0, 13, 4, 17};
	knownDFill[15] = {15, 0, 6};
	knownDFill[16] = {16, 0, 1, 4, 17, 5, 18};
	knownDFill[17] = {17, 0, 4, 14, 1, 16};
	knownDFill[18] = {18, 0, 1, 16, 5, 6};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 105);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 8.942685e+107, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 2.2618, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-4.47134e+107, -0.00379371, -0.00410766,
		-0.00246059, 1.11784e+107, -0.00318435, 0.158612, 0.0529346,
		-0.00516169, -0.00346739, -0.00302979, -0.00270841, -0.00226261,
		-0.00210017, -4.47134e+107, 0.110484, -0.118683, -1.3414e+107,
		-0.125485, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-2.61554, -0.409267, 0, -0.0738305,
		-0.0855092, 8.94269e+107, -1.3414e+108, -0.430117, 0, -0.110482,
		-0.064144, -0.069715, -0.0773028, -0.0803359, -0.0830468, -0.0877732,
		-0.0898745, -0.0918394, 0.000168841, -0.0256749, -0.00379539,
		-0.0383688, 0.0241089, 0.0361633, 4.44896e-05, -6.67345e-05, -0.0744708,
		-0.00408481, 0.0700659, 0.000297336, -0.0860146, -0.00245026, 0.0833116,
		0.000242396, -2.23567e+107, -0.0372559, 3.35351e+107, 0.0230163,
		4.78367e-05, -7.1755e-05, -0.0269713, -0.0032494, 0.025627, 0.0384405,
		-0.0403234, 4.28086e-05, -0.110824, 0.158615, 0.107707, 0.161561,
		0.000168841, -0.0653302, 0.170893, 0.000350842, -0.0704167, -0.00513103,
		0.0649348, 0.000320188, -0.0778975, -0.00344947, 0.0741507, 0.000279413,
		-0.0808948, -0.00301523, 0.0776001, 0.000264854, -0.0835765,
		-0.00269627, 0.0806154, 0.000252718, -0.088258, -0.00225369, 0.0857619,
		0.000233478, -0.0903414, -0.00209236, 0.0880156, 0.000225672,
		-0.0922907, -0.00195713, 0.0901079, 8.94269e+107, -1.3414e+108,
		-0.000168841, 0.110653, 0.110653, -0.0770246, -0.118455, 0.030695,
		-0.0289306, -0.043396, 5.33876e-05, -8.00813e-05, -4.02421e+107,
		-0.0710155, 2.68281e+107, -0.0276195, 5.7404e-05, -8.6106e-05,
		-0.0809138, -0.125349, -0.0307524, -0.0461286, 0.0322587, 5.13703e-05};
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
				XOLOTL_REQUIRE_CLOSE(hPartials[startingIdx + index],
					knownPartials[startingIdx + j], 0.01);
			}
			startingIdx += row.size();
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
