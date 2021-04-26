#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the Alloy network.
 */
BOOST_AUTO_TEST_SUITE(AlloyNetwork_testSuite)

BOOST_AUTO_TEST_CASE(fullyRefined)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 6 4" << std::endl
			  << "process=reaction sink" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = AlloyReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxSize = opts.getMaxImpurity();
	NetworkType network(
		{maxV, maxSize, maxSize, maxI, maxSize, maxSize}, 1, opts);

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 23);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 23);
	// TODO: check it is within a given range?
	auto deviceMemorySize = network.getDeviceMemorySize();
	BOOST_REQUIRE(deviceMemorySize > 120000);
	BOOST_REQUIRE(deviceMemorySize < 140000);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.011664, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 23);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 6);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 6);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 1);

	// Check the single vacancy
	auto vacancy = network.getSingleVacancy();
	BOOST_REQUIRE_EQUAL(vacancy.getId(), 17);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 9, 10, 11, 17, 18, 19, 20, 1, 2, 3, 12};
	knownDFill[1] = {1, 0, 9, 10, 17, 18, 19, 20, 21, 2, 3, 12, 11};
	knownDFill[2] = {2, 0, 10, 1, 9, 17, 18, 19, 20, 21, 22, 3, 11, 12};
	knownDFill[3] = {3, 0, 11, 1, 10, 2, 9, 17, 18, 19, 20, 21, 22, 4, 12};
	knownDFill[4] = {
		4, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 5, 6, 7, 8, 20};
	knownDFill[5] = {5, 4, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 6, 7, 8};
	knownDFill[6] = {6, 4, 10, 5, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 7, 8};
	knownDFill[7] = {
		7, 4, 11, 5, 10, 6, 9, 13, 14, 16, 17, 18, 19, 20, 21, 22, 8};
	knownDFill[8] = {
		8, 4, 12, 5, 11, 6, 10, 7, 9, 13, 14, 17, 18, 19, 20, 21, 22};
	knownDFill[9] = {9, 0, 20, 1, 21, 2, 22, 4, 19, 5, 6, 7, 8, 13, 10, 11, 12,
		14, 15, 16, 17, 18};
	knownDFill[10] = {10, 0, 19, 1, 20, 2, 21, 3, 22, 4, 18, 5, 6, 7, 8, 9, 11,
		12, 13, 14, 15, 16, 17};
	knownDFill[11] = {11, 0, 18, 1, 19, 2, 20, 3, 21, 4, 17, 5, 6, 7, 8, 9, 10,
		12, 13, 14, 15, 16, 22};
	knownDFill[12] = {12, 0, 17, 1, 18, 2, 19, 3, 20, 4, 5, 6, 7, 8, 9, 11, 10,
		13, 14, 15, 16, 21, 22};
	knownDFill[13] = {13, 4, 5, 6, 7, 8, 9, 14, 10, 11, 12, 17};
	knownDFill[14] = {14, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17};
	knownDFill[15] = {15, 4, 5, 6, 9, 16, 10, 11, 12, 17, 22, 18, 21, 19, 20};
	knownDFill[16] = {
		16, 4, 5, 6, 7, 9, 10, 11, 12, 15, 17, 18, 22, 19, 21, 20};
	knownDFill[17] = {17, 0, 1, 2, 3, 4, 5, 6, 13, 15, 7, 14, 16, 8, 9, 18, 10,
		19, 11, 20, 12, 21, 22};
	knownDFill[18] = {18, 0, 1, 2, 3, 4, 5, 13, 15, 6, 14, 16, 7, 8, 9, 19, 10,
		20, 11, 21, 12, 22, 17};
	knownDFill[19] = {19, 0, 1, 2, 3, 4, 13, 15, 5, 14, 16, 6, 7, 8, 9, 20, 10,
		21, 11, 22, 12, 17, 18};
	knownDFill[20] = {20, 0, 1, 2, 3, 4, 14, 16, 5, 6, 7, 8, 9, 21, 10, 22, 11,
		13, 15, 12, 17, 19, 18};
	knownDFill[21] = {21, 1, 2, 3, 6, 7, 8, 9, 22, 10, 13, 15, 11, 14, 16, 12,
		17, 20, 18, 19};
	knownDFill[22] = {
		22, 2, 3, 7, 8, 9, 13, 15, 10, 14, 16, 11, 12, 17, 21, 18, 20, 19};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 413);
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
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 5048430879.96, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::Void, 1), 2.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::Perfect, 3), 30.0,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Frank, 6), 1.394446,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Faulted, 1), 0.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 7782153079, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-1.79219e+09, 1.85298e+09, 3.81785e+09,
		9.93479e+09, -9.45535e+09, -4.17172e+09, -7.68927e+08, 2.88205e+09,
		8.58606e+09, -7.85155e+10, -3.45865e+10, -2.27619e+10, -1.39981e+10,
		-5.3704e+09, -1.02448e+10, -3.69144e+09, -6.92805e+09, 4.06728e+08,
		2.84616e+08, 1.50115e+09, 3.00266e+09, 3.1405e+09, 4.93097e+09, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-7.78215e+09, -5.78851e+08,
		2.67961e+08, 9.4625e+08, 12438.8, 11978.4, 11558.6, -87954.2, 333115,
		179225, 126275, 3.56182e+09, -6.54255e+09, 4.14067e+09, -1.60614e+08,
		4.55379e+06, 11321.8, 10978.6, -118915, -91113.9, -74190.8, 344437,
		184493, 2.24503e+09, 4.01861e+09, -4.44836e+09, 2.15954e+09,
		2.15954e+09, 4.30128e+09, -1.4619e+08, 10435.7, -179225, -122742,
		-93995.5, -76503, -64702.4, 354872, 1.80502e+09, 1.80502e+09, -907422,
		1.48125e+09, 1.48125e+09, 2.24047e+09, 2.24047e+09, 4.44747e+09,
		4.44747e+09, -354872, -184493, -126275, -96656.3, -78638.7, -66487,
		3.58296e+08, 3.17469e+09, -1.03899e+10, -3.76754e+09, -2.02102e+09,
		-1.42287e+09, -1.11839e+09, -1.05684e+08, -1.08527e+08, -8.78048e+07,
		-8.90013e+07, -1.25903e+07, -2.23786e+07, -3.00636e+07, 6.29138e+07,
		5.65822e+07, 5.13656e+07, 4.71167e+07, 4.71167e+07, -8.08967e+09,
		3.76754e+09, -1.34104e+08, -2.08005e+09, -1.45647e+09, -8.72402e+07,
		-8.9526e+07, -7.27232e+07, -7.36831e+07, -7.15151e+07, -1.56958e+07,
		-2.1501e+07, -6.9211e+07, 5.43125e+07, 4.99631e+07, 4.61292e+07,
		-6.78453e+09, 2.02102e+09, -1.15428e+08, 3.90164e+09, -1.23463e+08,
		-7.47239e+07, -7.66364e+07, -6.24591e+07, -6.32614e+07, -6.27866e+06,
		-6.82595e+07, -5.82218e+07, -5.9539e+07, -6.06558e+07, 4.80338e+07,
		4.49048e+07, -4.68922e+09, 1.42287e+09, 1.42287e+09, 2.08005e+09,
		2.08005e+09, 4.02511e+09, -1.14538e+08, -6.56409e+07, -6.72857e+07,
		-5.56809e+07, -4.80237e+06, -4.99631e+07, -1.02731e+08, -5.24944e+07,
		-5.34523e+07, -5.42917e+07, 4.32314e+07, -4.44048e+08, 1.47669e+09,
		1.11839e+09, 1.45647e+09, 1.45647e+09, 2.13645e+09, 2.13645e+09,
		4.13965e+09, 4.13965e+09, -5.87301e+07, -6.01734e+07, -4.32314e+07,
		-4.49048e+07, -4.61292e+07, -9.42333e+07, -4.79553e+07, -4.86906e+07,
		-9.2016e+10, -4.14058e+09, -2.78022e+09, -4.30121e+09, -2.85665e+09,
		-4.44741e+09, -2.92234e+09, -3.68611e+09, -1.91407e+09, -3.83243e+09,
		-3.96445e+09, -4.08535e+09, 5.87301e+07, -4.38874e+09, -2.66788e+09,
		-2.72098e+09, -2.78814e+09, -4.58222e+09, -3.0297e+09, -3.07791e+09,
		-1.19014e+09, -1.69712e+09, -4.72394e+10, -2.15943e+09, -1.38885e+09,
		-2.24038e+09, -1.43713e+09, 76503, -1.47715e+09, 66487, -1.51157e+09,
		-1.94206e+09, -5.78712e+08, -2.01242e+09, -2.07691e+09, 5.34523e+07,
		4.86906e+07, 6.09466e+09, -1.50422e+09, -1.49655e+09, -2.31421e+09,
		-2.38223e+09, -1.58681e+09, -1.61091e+09, -4.11028e+08, -3.17956e+10,
		-1.48108e+09, -9.0531e+08, 118915, -9.46405e+08, 93995.5, -9.78925e+08,
		78638.7, -1.00609e+09, -1.34737e+09, -1.35153e+08, -1.39081e+09,
		5.82218e+07, 5.24944e+07, 4.79553e+07, 3.10083e+08, 1.5747e+09,
		-1.09239e+09, -1.58488e+09, -1.6305e+09, -1.09152e+09, -1.10759e+09,
		-1.07389e+09, -2.38957e+10, 320676, -6.49394e+08, 173515, -6.91726e+08,
		122742, -7.22189e+08, 96656.3, -7.46534e+08, -1.11839e+09, 6.29138e+07,
		5.65822e+07, 5.13656e+07, 4.71167e+07, 1.30296e+08, 1.88709e+09,
		3.37907e+09, -1.2137e+09, -1.24806e+09, -8.38734e+08, -8.50787e+08,
		-8.10701e+08, -8.25523e+08, -9.95263e+09, -1.05684e+08, -8.72402e+07,
		-7.47239e+07, -6.56409e+07, -5.87301e+07, 1.34749e+08, 4.58223e+09,
		-2.31421e+09, -1.58488e+09, -1.2137e+09, -344437, -1.02452e+10,
		-1.08527e+08, -8.9526e+07, -7.66364e+07, -6.72857e+07, -6.01734e+07,
		-4.58222e+09, -2.38223e+09, -1.6305e+09, -1.24806e+09, 344437, 344437,
		-6.77003e+09, -8.78048e+07, -7.27232e+07, -6.24591e+07, 4.82093e+07,
		3.07791e+09, -1.58681e+09, -1.09152e+09, -8.38734e+08, 41175.4, 322739,
		199146, 199146, 167749, 167749, -6.92883e+09, -8.90013e+07,
		-7.36831e+07, -6.32614e+07, -5.56809e+07, -3.07791e+09, -1.61091e+09,
		-1.10759e+09, -8.50787e+08, 281564, 281564, 193339, 193339, 156734,
		156734, 295039, -5.84714e+09, -320676, -333115, -344437, -354872,
		-7.55041e+07, -1.25828e+08, 8.28706e+07, 7.43795e+07, 6.21804e+07,
		7.49328e+07, 6.72882e+07, 5.56823e+07, -4.32314e+07, 1.43806e+08,
		2.69309e+09, 1.22277e+08, 1.45667e+09, 1.08228e+08, 1.03124e+09,
		9.80725e+07, 8.10391e+08, -316794, -6.16223e+09, -167247, -173515,
		-179225, -184493, -7.89608e+07, 9.43045e+07, 8.72402e+07, 7.27232e+07,
		2.67334e+07, 7.66364e+07, 6.32614e+07, -4.99631e+07, -4.49048e+07,
		1.00916e+08, 2.76937e+09, 9.06005e+07, 1.49655e+09, 8.29874e+07,
		1.05392e+09, 7.70413e+07, 8.2533e+08, 566405, -6.36287e+09, -114716,
		-118915, -122742, -126275, 1.12059e+08, 1.31938e+09, 9.26539e+08,
		9.55788e+07, 8.9526e+07, 7.36831e+07, -5.82218e+07, -1.02731e+08,
		-4.61292e+07, 8.03497e+07, 2.84941e+09, 7.40857e+07, 1.53052e+09,
		6.91443e+07, 1.07389e+09, 1.27876e+09, 28817.5, 146344, -6.44874e+09,
		-87954.2, -91113.9, -93995.5, -96656.3, 1.97528e+08, 1.35658e+09,
		9.39788e+08, -6.9211e+07, -5.9539e+07, -5.24944e+07, -9.42333e+07,
		6.78568e+07, 2.91739e+09, 6.35758e+07, 1.56033e+09, 1.64489e+09,
		1.58488e+09, 1.09152e+09, 1.3051e+09, 12109, 175376, 314032,
		-6.47587e+09, -74190.8, -76503, -78638.7, -6.06558e+07, -5.34523e+07,
		-4.79553e+07, 5.9317e+07, 2.9767e+09, 2.37034e+09, 2.31421e+09,
		1.58681e+09, 1.68397e+09, 1.6305e+09, 1.10759e+09, -8.10701e+08,
		5671.81, 331016, 26452.3, 68864.3, -6.54007e+09, -64702.4, -66487,
		-5.42917e+07, -4.86906e+07, 4.50048e+09, 4.44747e+09, 3.0297e+09,
		2.43282e+09, 2.38223e+09, 1.61091e+09, -1.07389e+09, -8.25523e+08,
		2605.09, 325344, 15301.8, 208641, 373467};
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
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 17);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.14068, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 1.9, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 29106.01, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 103680000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto& reg0 = cluster.getRegion();
	Composition lo = reg0.getOrigin();
	Composition hi = reg0.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::V] = 0;
	comp[Spec::I] = 2;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 10);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.17725, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 5.4559, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 156594191.2, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 51840000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.5, 0.01);
	const auto& reg1 = cluster.getRegion();
	lo = reg1.getOrigin();
	hi = reg1.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::I] = 0;
	comp[Spec::Void] = 7;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 15);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.26912, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 10.9416, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg2 = cluster.getRegion();
	lo = reg2.getOrigin();
	hi = reg2.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 7);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 8);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Void] = 0;
	comp[Spec::Frank] = 8;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 3);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.49776, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 9.4, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg3 = cluster.getRegion();
	lo = reg3.getOrigin();
	hi = reg3.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 8);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 9);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Frank] = 0;
	comp[Spec::Perfect] = 5;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 5);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.321142, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.2480, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 6150067.77, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 20736000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.7, 0.01);
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 6);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Perfect] = 0;
	comp[Spec::Faulted] = 7;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 13);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.4656, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.2186, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg5 = cluster.getRegion();
	lo = reg5.getOrigin();
	hi = reg5.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 7);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 8);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
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
	paramFile << "netParam=20 0 0 6 4" << std::endl
			  << "grouping=10 5" << std::endl
			  << "process=reaction sink" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = AlloyReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxSize = opts.getMaxImpurity();
	NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
	// Adapt maxSize
	int i = 0;
	while (maxSize + 1 > pow(groupingWidth, i)) {
		++i;
	}
	maxSize = pow(groupingWidth, i) - 1;

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxV, maxSize, maxSize, maxI, maxSize, maxSize};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {{maxV + 1,
		groupingWidth, groupingWidth, maxI + 1, groupingWidth, groupingWidth}};
	NetworkType network = NetworkType(maxSpeciesAmounts, subdivRatios, 1, opts);

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 47);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 53);
	// TODO: check it is within a given range?
	auto deviceMemorySize = network.getDeviceMemorySize();
	BOOST_REQUIRE(deviceMemorySize > 700000);
	BOOST_REQUIRE(deviceMemorySize < 1100000);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.011664, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 47);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 6);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 6);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 1);

	// Check the single vacancy
	auto vacancy = network.getSingleVacancy();
	BOOST_REQUIRE_EQUAL(vacancy.getId(), 21);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 2, 47, 3, 48, 5, 6, 7, 9, 10, 11, 12, 13, 14, 49, 15,
		50, 16, 17, 18, 51, 19, 52, 20, 21, 22, 23, 27, 28, 29, 30, 31, 32, 33,
		34, 35, 36, 37, 38, 39, 43, 44, 45, 46, 26, 8, 40, 41, 42, 24, 25};
	knownDFill[1] = {1, 0, 5, 6, 7, 9, 10, 11, 12, 21, 22, 23, 24, 31, 32, 33,
		34, 35, 36, 37, 38, 2, 47, 25, 26, 27, 28, 29, 30};
	knownDFill[2] = {2, 0, 47, 6, 27, 28, 29, 30, 31, 32, 33, 34, 1, 5, 9, 10,
		11, 12, 21, 22, 23, 24, 25, 26, 3, 48, 4};
	knownDFill[3] = {3, 0, 2, 47, 48, 7, 35, 36, 37, 38, 1, 6, 5, 9, 10, 11, 12,
		31, 32, 33, 34, 21, 22, 23, 24, 25, 26, 4, 27, 28, 29, 30};
	knownDFill[4] = {4, 0, 3, 48, 39, 1, 7, 2, 47, 6, 31, 32, 33, 34, 5, 9, 10,
		11, 12, 21, 22, 23, 24, 25, 26, 27, 38, 28, 37, 29, 36, 30, 35};
	knownDFill[5] = {5, 0, 9, 1, 2, 47, 3, 48, 6, 7, 10, 11, 12, 13, 14, 49, 15,
		50, 16, 17, 18, 51, 19, 52, 20, 21, 22, 23, 24, 27, 28, 29, 30, 31, 32,
		33, 34, 35, 36, 37, 38, 43, 44, 45, 46, 25, 8, 39, 40, 41, 42, 26};
	knownDFill[6] = {6, 0, 31, 1, 2, 47, 5, 9, 10, 11, 12, 13, 14, 49, 15, 50,
		16, 18, 51, 19, 52, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33,
		34, 43, 44, 7, 8, 40, 35, 36, 37, 38, 39, 41, 42};
	knownDFill[7] = {7, 0, 35, 1, 5, 6, 9, 10, 11, 12, 13, 14, 49, 15, 50, 16,
		19, 52, 20, 21, 22, 23, 24, 25, 26, 43, 44, 8, 38, 37, 36, 39, 40, 41,
		42, 31, 34, 32, 33};
	knownDFill[8] = {8, 0, 39, 5, 7, 6, 13, 14, 49, 15, 50, 16, 21, 22, 23, 24,
		25, 26, 43, 44, 9, 42, 10, 41, 11, 40, 12, 31, 38, 32, 37, 33, 36, 34,
		35};
	knownDFill[9] = {9, 0, 23, 1, 24, 2, 47, 3, 48, 5, 6, 44, 7, 14, 49, 8, 15,
		50, 10, 11, 12, 13, 16, 17, 18, 51, 19, 52, 20, 21, 22, 25, 26, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46};
	knownDFill[10] = {10, 0, 22, 1, 23, 2, 47, 3, 48, 5, 6, 43, 7, 14, 49, 8,
		15, 50, 9, 11, 12, 13, 16, 17, 18, 51, 19, 52, 20, 21, 24, 25, 26, 27,
		28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 44, 45, 46, 42};
	knownDFill[11] = {11, 0, 21, 1, 22, 2, 47, 3, 48, 5, 6, 13, 7, 14, 49, 8,
		15, 50, 9, 10, 12, 16, 17, 18, 51, 19, 52, 20, 23, 24, 25, 26, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 43, 44, 45, 46, 41, 42};
	knownDFill[12] = {12, 0, 1, 21, 2, 47, 26, 3, 48, 5, 6, 7, 14, 49, 8, 15,
		50, 9, 11, 10, 13, 16, 17, 18, 51, 19, 52, 20, 22, 23, 24, 25, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 43, 44, 45, 46, 40, 41, 42};
	knownDFill[13] = {13, 0, 14, 49, 5, 6, 15, 50, 7, 8, 9, 43, 10, 44, 11, 12,
		21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
		16};
	knownDFill[14] = {14, 0, 49, 15, 50, 5, 6, 16, 7, 8, 9, 10, 11, 12, 13, 23,
		24, 25, 26, 21, 22, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44,
		43};
	knownDFill[15] = {15, 0, 50, 16, 5, 6, 7, 8, 9, 10, 11, 12, 14, 49, 21, 22,
		23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44};
	knownDFill[16] = {16, 0, 5, 6, 7, 8, 9, 10, 11, 12, 14, 49, 26, 15, 50, 21,
		22, 23, 24, 25, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42};
	knownDFill[17] = {17, 0, 18, 51, 5, 6, 19, 52, 9, 45, 10, 46, 11, 12, 21,
		22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 20, 37};
	knownDFill[18] = {18, 0, 51, 19, 52, 5, 6, 20, 9, 10, 11, 12, 17, 23, 24,
		25, 26, 21, 22, 31, 32, 33, 34, 35, 36, 37, 46, 45};
	knownDFill[19] = {19, 0, 52, 20, 5, 6, 7, 9, 10, 11, 12, 18, 51, 21, 22, 23,
		24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 46};
	knownDFill[20] = {20, 0, 5, 6, 7, 9, 10, 11, 12, 18, 51, 26, 19, 52, 21, 22,
		23, 24, 25, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42};
	knownDFill[21] = {21, 0, 1, 2, 47, 3, 48, 4, 5, 6, 14, 49, 18, 51, 7, 15,
		50, 19, 52, 8, 9, 22, 10, 23, 11, 24, 12, 25, 13, 31, 34, 35, 36, 37,
		38, 39, 40, 41, 16, 42, 17, 20, 26, 27, 28, 29, 30, 32, 33, 43, 44, 45,
		46};
	knownDFill[22] = {22, 0, 1, 2, 47, 3, 48, 4, 5, 13, 17, 6, 14, 49, 18, 51,
		7, 15, 50, 19, 52, 8, 9, 23, 10, 24, 11, 25, 12, 26, 33, 34, 35, 36, 37,
		38, 39, 40, 16, 41, 20, 21, 27, 28, 29, 30, 31, 32, 42, 43, 44, 45, 46};
	knownDFill[23] = {23, 0, 13, 17, 1, 2, 47, 3, 48, 4, 5, 43, 45, 6, 14, 49,
		18, 51, 7, 15, 50, 19, 52, 8, 9, 24, 10, 25, 11, 26, 12, 32, 33, 34, 35,
		36, 37, 38, 39, 16, 40, 20, 21, 22, 27, 28, 29, 30, 31, 41, 42, 44, 46};
	knownDFill[24] = {24, 0, 43, 45, 1, 2, 47, 3, 48, 4, 5, 44, 46, 6, 14, 49,
		18, 51, 7, 15, 50, 19, 52, 8, 9, 25, 10, 26, 11, 13, 17, 12, 31, 32, 33,
		34, 35, 36, 37, 38, 16, 39, 20, 21, 23, 22, 27, 28, 29, 30, 40, 41, 42};
	knownDFill[25] = {25, 0, 44, 46, 2, 47, 3, 48, 4, 5, 14, 49, 18, 51, 6, 15,
		50, 19, 52, 7, 16, 20, 8, 9, 26, 10, 13, 17, 11, 43, 45, 12, 31, 32, 33,
		34, 35, 36, 37, 38, 21, 24, 22, 23, 27, 28, 29, 30, 39, 40, 41, 42};
	knownDFill[26] = {26, 0, 14, 49, 18, 51, 2, 47, 3, 48, 4, 5, 6, 15, 50, 19,
		52, 7, 8, 9, 13, 17, 10, 43, 45, 11, 44, 46, 12, 31, 32, 33, 34, 35, 36,
		37, 16, 38, 20, 21, 25, 22, 24, 23, 28, 29, 30, 39, 40, 41, 42};
	knownDFill[27] = {27, 0, 1, 9, 2, 47, 24, 25, 26, 5, 6, 10, 12, 11, 21, 28,
		22, 29, 23, 30, 31, 32, 33, 34, 35, 36, 37, 38};
	knownDFill[28] = {28, 0, 1, 10, 2, 47, 23, 24, 25, 26, 5, 6, 9, 27, 11, 12,
		21, 29, 22, 30, 31, 32, 33, 34, 35, 36, 37};
	knownDFill[29] = {29, 0, 1, 11, 2, 47, 22, 23, 24, 25, 26, 5, 6, 9, 28, 10,
		27, 12, 21, 30, 31, 32, 33, 34, 35, 36};
	knownDFill[30] = {30, 0, 1, 5, 12, 2, 47, 21, 22, 23, 24, 25, 3, 48, 26, 6,
		9, 29, 10, 28, 11, 27, 31, 32, 33, 34, 35};
	knownDFill[31] = {31, 0, 10, 1, 2, 47, 5, 9, 6, 24, 7, 44, 8, 14, 49, 11,
		12, 13, 37, 39, 40, 41, 42, 15, 50, 16, 17, 18, 51, 19, 52, 20, 21, 32,
		22, 33, 23, 34, 25, 35, 26, 36, 27, 28, 29, 30, 38, 43, 45, 46};
	knownDFill[32] = {32, 0, 11, 1, 2, 47, 5, 10, 6, 23, 7, 43, 8, 14, 49, 9,
		31, 12, 13, 38, 40, 41, 42, 15, 50, 16, 18, 51, 19, 52, 20, 21, 33, 22,
		34, 24, 35, 25, 36, 26, 37, 27, 28, 29, 30, 44, 45, 46, 39};
	knownDFill[33] = {33, 0, 12, 1, 2, 47, 5, 11, 6, 22, 7, 13, 8, 14, 49, 9,
		32, 10, 31, 41, 42, 15, 50, 16, 18, 51, 19, 52, 20, 21, 34, 23, 35, 24,
		36, 25, 37, 26, 38, 27, 28, 29, 30, 43, 44, 46, 39, 40};
	knownDFill[34] = {34, 0, 5, 1, 2, 47, 12, 6, 21, 7, 26, 8, 14, 49, 9, 33,
		10, 32, 11, 31, 13, 39, 42, 15, 50, 16, 18, 51, 19, 52, 20, 22, 35, 23,
		36, 24, 37, 25, 38, 27, 28, 29, 30, 43, 44, 40, 41};
	knownDFill[35] = {35, 0, 32, 1, 5, 31, 6, 9, 7, 24, 8, 44, 10, 34, 11, 33,
		12, 13, 41, 14, 49, 15, 50, 16, 18, 51, 19, 52, 20, 21, 36, 22, 37, 23,
		38, 25, 39, 26, 40, 27, 28, 29, 30, 43, 42};
	knownDFill[36] = {36, 0, 33, 1, 5, 32, 6, 10, 7, 23, 8, 43, 9, 35, 11, 34,
		12, 13, 42, 14, 49, 15, 50, 16, 18, 51, 19, 52, 20, 21, 37, 22, 38, 24,
		39, 25, 40, 26, 41, 27, 28, 29, 31, 44};
	knownDFill[37] = {37, 0, 34, 1, 5, 33, 6, 11, 7, 22, 8, 13, 9, 36, 10, 35,
		12, 14, 49, 15, 50, 16, 18, 51, 19, 52, 20, 21, 38, 23, 39, 24, 40, 25,
		41, 26, 42, 27, 28, 31, 32, 43, 44};
	knownDFill[38] = {38, 0, 6, 1, 5, 34, 12, 7, 21, 8, 26, 9, 37, 10, 36, 11,
		35, 13, 14, 49, 15, 50, 16, 19, 52, 20, 22, 39, 23, 40, 24, 41, 25, 42,
		27, 31, 33, 32, 43, 44};
	knownDFill[39] = {39, 0, 36, 5, 35, 6, 31, 7, 9, 8, 24, 10, 38, 11, 37, 12,
		13, 14, 49, 15, 50, 16, 19, 52, 20, 21, 40, 22, 41, 23, 42, 25, 26, 32,
		34, 33, 43, 44};
	knownDFill[40] = {40, 0, 37, 5, 36, 6, 32, 7, 10, 8, 23, 9, 39, 11, 38, 12,
		13, 14, 49, 15, 50, 16, 19, 52, 20, 21, 41, 22, 42, 24, 25, 26, 31, 35,
		33, 34, 43, 44};
	knownDFill[41] = {41, 0, 38, 5, 37, 6, 33, 7, 11, 8, 22, 9, 40, 10, 39, 12,
		13, 14, 49, 15, 50, 16, 19, 52, 20, 21, 42, 23, 24, 25, 26, 31, 36, 32,
		35, 34, 43, 44};
	knownDFill[42] = {42, 0, 7, 5, 38, 6, 34, 12, 8, 21, 9, 41, 10, 40, 11, 39,
		13, 14, 49, 15, 50, 16, 20, 22, 23, 24, 25, 26, 31, 37, 32, 36, 33, 35,
		43, 44};
	knownDFill[43] = {43, 0, 14, 49, 5, 6, 15, 50, 7, 8, 9, 44, 10, 11, 12, 13,
		21, 31, 32, 33, 34, 35, 16, 36, 22, 23, 24, 25, 26, 37, 38, 39, 40, 41,
		42};
	knownDFill[44] = {44, 0, 14, 49, 5, 6, 15, 50, 7, 8, 9, 10, 11, 12, 13, 22,
		31, 32, 33, 34, 16, 35, 21, 43, 23, 24, 25, 26, 36, 37, 38, 39, 40, 41,
		42};
	knownDFill[45] = {45, 0, 18, 51, 5, 6, 19, 52, 9, 46, 10, 11, 12, 17, 21,
		31, 32, 33, 34, 35, 20, 36, 22, 26, 23, 25, 24};
	knownDFill[46] = {46, 0, 18, 51, 5, 6, 19, 52, 9, 10, 11, 12, 17, 22, 31,
		32, 33, 34, 20, 35, 21, 45, 23, 26, 24, 25};
	knownDFill[47] = {47, 0, 2, 6, 27, 28, 29, 30, 31, 32, 33, 34, 1, 5, 9, 10,
		11, 12, 21, 22, 23, 24, 25, 26, 3, 48, 4};
	knownDFill[48] = {48, 0, 2, 47, 3, 7, 35, 36, 37, 38, 1, 6, 5, 9, 10, 11,
		12, 31, 32, 33, 34, 21, 22, 23, 24, 25, 26, 4, 27, 28, 29, 30};
	knownDFill[49] = {49, 0, 14, 15, 50, 5, 6, 16, 7, 8, 9, 10, 11, 12, 13, 23,
		24, 25, 26, 21, 22, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44,
		43};
	knownDFill[50] = {50, 0, 15, 16, 5, 6, 7, 8, 9, 10, 11, 12, 14, 49, 21, 22,
		23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44};
	knownDFill[51] = {51, 0, 18, 19, 52, 5, 6, 20, 9, 10, 11, 12, 17, 23, 24,
		25, 26, 21, 22, 31, 32, 33, 34, 35, 36, 37, 46, 45};
	knownDFill[52] = {52, 0, 19, 20, 5, 6, 7, 9, 10, 11, 12, 18, 51, 21, 22, 23,
		24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 46};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 2055);
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
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 5735953497.9, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::Void, 1), 14.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::Perfect, 3), 204.0,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Frank, 6), 9.377118,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Faulted, 1), 0.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 16881508482, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-1.60522e+10, -3.5809e+09, -5.17624e+09,
		5.73588e+09, 2.54063e+10, -1.1317e+10, -3.25617e+09, -1.26326e+09,
		1.1073e+10, -2.99435e+11, -1.42494e+11, -9.38683e+10, -6.78426e+10,
		-4.43316e+09, -6.37309e+08, 2.46714e+09, -1.33153e+10, -2.61011e+09,
		-1.14603e+08, 1.52464e+09, -8.23282e+09, 9.28303e+08, 8.57322e+08,
		2.1323e+09, 3.70545e+09, 6.03226e+09, 6.0872e+09, -1.49871e+09,
		-1.91548e+09, -5.03716e+08, -8.39181e+08, -8.65364e+09, -6.266e+09,
		-4.22522e+09, -3.69568e+09, -3.07372e+09, -2.64713e+09, -2.25522e+09,
		-1.72707e+09, -9.63242e+08, 7.1752e+08, 2.76546e+09, 5.50675e+09,
		-6.41722e+09, -1.06913e+10, -3.92022e+09, -6.6945e+09, -1.03656e+10,
		-9.11985e+09, -1.40574e+10, -1.39882e+10, -8.90744e+09, -8.73621e+09,
		0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-1.75722e+10, -9.92391e+07,
		-5.90879e+08, 0, -2.5611e+08, 3.84165e+08, -2.70728e+08, -2.52393e+08,
		-2.1052e+08, -3.76754e+09, -2.02102e+09, -1.42287e+09, -1.11839e+09,
		-6.04921e+07, -4.06145e+08, 1.87448e+07, -5.73286e+08, 1.01175e+08,
		-1.3321e+08, -8.78048e+07, -4.63596e+08, 0, -4.81724e+08, 0,
		-9.81908e+07, -1.25903e+07, -2.23786e+07, -3.00636e+07, -1.02608e+08,
		-1.05684e+08, -1.08527e+08, -1.11179e+08, -2.61644e+08, -2.56695e+08,
		-2.54027e+08, -2.52756e+08, -2.45628e+08, -2.46229e+08, -2.472e+08,
		-1.65235e+08, -2.51354e+08, -6.54094e+07, -6.99251e+07, -8.90013e+07,
		-9.00699e+07, 4.06739e+07, 6.56029e+07, 0, -3.34101e+07, 1.02561e+08,
		4.71167e+07, 4.36113e+07, -9.57088e+09, -9.92391e+07, -8.20626e+07,
		-4.60708e+07, -3.31446e+07, -5.78851e+08, 2.67961e+08, 9.4625e+08,
		2.42613e+09, 12438.8, 11978.4, 11558.6, 11182.3, -7.03943e+07,
		-6.19191e+07, -5.54202e+07, -5.02548e+07, -4.26074e+07, -3.96893e+07,
		-3.71943e+07, -3.50344e+07, 158413, -244289, -85875.6, 0, 333115,
		179225, 126275, 99136.5, -4.62179e+09, 9.5475e+07, -5.12029e+09,
		4.01038e+07, 3.11512e+08, 6.26374e+08, 1.1071e+09, 2.06182e+09,
		1.41585e+08, 1.42552e+08, 1.04507e+08, 7.78455e+07, 6.40102e+07,
		1.61312e+08, -2.0805e+09, -1.64654e+09, -1.16412e+09, -6.5361e+08,
		-7135.3, -3619.67, 0, 3675.96, 7389.5, 38316.8, 511674, -407229,
		16388.8, -2.40084e+09, 2.86265e+08, 4.55449e+09, 5.10625e+09,
		-2.40334e+09, 2.87909e+07, 1.55701e+08, 1.23121e+08, 9.18856e+07,
		6.08667e+07, 4.01192e+07, 1.61496e+08, 2.45475e+08, -2.764e+08,
		1.47378e+09, 1.75339e+09, 1.63728e+09, 1.65877e+08, 1.38475e+08,
		1.54375e+08, 1.61407e+08, 173994, 90036.2, 31345.7, -22014.9, -73842.5,
		-78721.4, 209498, 4.4533e+07, 4.9628e+07, 5.5668e+07, 6.31309e+07,
		-1.12943e+06, 1.44041e+08, 1.20015e+10, 1.20189e+10, 1.44041e+08,
		3.31446e+07, 1.58006e+08, 3.34046e+08, 7.22383e+07, 1.39846e+08,
		3.62243e+08, 2.50853e+08, 1.65358e+08, 9.76601e+07, 1.96054e+07,
		1.64921e+10, 5.69158e+09, 1.94202e+09, 0, -444223, -229787, -156728,
		-119638, -97113.2, -81943.9, 3.60995e+07, 1.48924e+08, 3.93892e+07,
		1.44478e+08, 4.31172e+07, 1.43388e+08, 4.72665e+07, 1.44927e+08,
		-1.63822e+10, 3.4339e+09, -1.34104e+08, -8.20626e+07, -4.86459e+08, 0,
		-1.05256e+08, 2.10512e+08, -2.19184e+08, -1.3364e+08, -2.08005e+09,
		-1.45647e+09, -1.13906e+09, -4.51776e+07, -3.10682e+08, 1.67718e+07,
		-4.9348e+08, 6.56029e+07, -1.09418e+08, -7.27232e+07, -3.84145e+08,
		2.98023e-08, -3.98694e+08, 0, -8.12207e+07, -7.15151e+07, -1.56958e+07,
		-2.1501e+07, -2.63456e+07, -8.47686e+07, -8.72402e+07, -8.9526e+07,
		-9.16596e+07, -2.35482e+08, -2.28035e+08, -2.23555e+08, -2.20811e+08,
		-2.18316e+08, -2.11544e+08, -2.11484e+08, -2.11755e+08, -4.9161e+07,
		-5.28245e+07, -7.36831e+07, -7.45656e+07, 4.00844e+07, 9.84043e+07,
		7.34543e+07, 3.4978e+07, 0, -3.19669e+07, 3.7695e+07, -1.35992e+10,
		-1.33953e+08, 1.0835e+09, -4.60708e+07, -5.3986e+07, 1.07972e+08,
		4.79164e+07, -9.50912e+07, -9.41276e+07, -9.16506e+07, -8.72705e+07,
		-1.71898e+07, -2.39271e+08, 6.13193e+07, -2.90308e+08, 0, -6.02003e+07,
		-1.72917e+08, -8.64583e+07, -2.24031e+08, 0, -4.56121e+07, -2.57051e+06,
		-4.98534e+06, -7.1743e+06, -9.15556e+06, -1.09546e+07, -5.32694e+07,
		-4.755e+07, -4.88725e+07, -5.00326e+07, -5.11169e+07, 1.30192e+09,
		2.03111e+09, 4.14655e+09, -1.90611e+07, -2.08315e+07, 2.91298e+07,
		9.1979e+07, 3.16828e+07, 3.37646e+07, 3.26633e+07, 3.14491e+07,
		3.02572e+07, 2.80783e+07, 3.09715e+07, 3.02854e+07, -1.20788e+10,
		-1.42499e+08, 1.41434e+09, -3.31446e+07, -1.20089e+08, 1.29634e+08,
		-7.58027e+07, -7.62098e+07, -7.59765e+07, -7.52277e+07, -3.50544e+07,
		-2.32166e+08, -3.86943e+07, -2.07228e+08, 0, -4.28928e+07, -1.29168e+08,
		-6.45841e+07, -3.28355e+07, -1.25763e+06, -2.4955e+06, -3.66449e+06,
		-4.76113e+06, -5.7886e+06, -2.95247e+07, -3.5903e+07, -3.66985e+07,
		2.33412e+07, 4.77977e+09, 2.46101e+09, 1.67016e+09, 2.54295e+07,
		2.50204e+07, 2.44919e+07, 2.39199e+07, 1.21073e+08, 1.21073e+08,
		1.17334e+08, 1.17334e+08, -6.38437e+08, 1.44041e+08, 1.50562e+09,
		1.24862e+08, 1.24862e+08, 1.93832e+08, -2.78165e+07, -1.53298e+08, 0,
		-1.96809e+08, -3.28014e+07, -3.38975e+07, -2.1559e+07, -2.21406e+07,
		-2.26125e+07, -2.30042e+07, -2.33412e+07, -2.36386e+07, -2.84585e+07,
		-2.90599e+07, 5.13327e+09, 5.13327e+09, 2.64217e+09, 2.64217e+09,
		1.7928e+09, 1.7928e+09, 1.36158e+09, 1.12825e+08, 1.12825e+08,
		1.05089e+08, 1.05089e+08, 1.0027e+08, 1.0027e+08, 9.76601e+07,
		9.76601e+07, -3.13484e+11, -3.68611e+09, -1.91407e+09, -4.14058e+09,
		-2.78022e+09, -2.51835e+10, 0, -2.74868e+10, 0, -3.83243e+09,
		-4.39126e+09, -4.65649e+09, -4.89427e+09, -2.49641e+10, 2.78758e+07,
		-3.28014e+07, -2.73074e+10, 1.62476e+07, -2.66788e+09, -2.72098e+09,
		-2.78814e+09, -4.38874e+09, -5.73595e+09, -3.0297e+09, -1.61809e+10,
		0.00610785, -1.69656e+10, 0.000922441, -3.47348e+09, -1.19014e+09,
		-1.69712e+09, -2.85665e+09, -2.92234e+09, -4.30121e+09, -4.44741e+09,
		-4.58222e+09, -4.70761e+09, -3.96445e+09, -4.08535e+09, -4.18794e+09,
		-4.29271e+09, -4.383e+09, -4.52582e+09, -4.65744e+09, -4.77977e+09,
		-4.81067e+09, -4.92426e+09, -5.03154e+09, -5.13327e+09, -4.52765e+09,
		-3.07791e+09, -3.12225e+09, -1.55746e+11, -1.94206e+09, -5.78712e+08,
		-2.15943e+09, -1.38885e+09, -1.30605e+10, 0, -1.13832e+10, 5.69158e+09,
		-2.01242e+09, -2.28478e+09, -2.3322e+09, -2.53722e+09, -1.28543e+10,
		2.41802e+07, 0, -1.40888e+10, 7.77443e+07, 6.09466e+09, -1.50422e+09,
		-1.49655e+09, -2.26093e+09, -2.96709e+09, -1.58681e+09, -8.45026e+09, 0,
		-8.8426e+09, 0, -1.8087e+09, -4.11028e+08, -1.43713e+09, -1.47715e+09,
		-1.51157e+09, -2.24038e+09, -2.31413e+09, -2.38217e+09, -2.44558e+09,
		-2.07691e+09, -2.13649e+09, -2.19199e+09, -2.2356e+09, -2.33141e+09,
		-2.2814e+09, -2.37406e+09, -2.45889e+09, -2.61008e+09, -2.49606e+09,
		-2.57074e+09, -2.39832e+09, -1.61091e+09, -1.63308e+09, 3.4187e+07,
		-1.0351e+11, -1.34737e+09, -1.35153e+08, -1.48108e+09, -9.0531e+08,
		-8.92415e+09, 0, -5.82607e+09, 5.82607e+09, -1.39081e+09, -1.56588e+09,
		-1.53601e+09, -1.70744e+09, -8.72935e+09, 2.1193e+07, 3.28014e+07,
		-9.60741e+09, 1.0562e+08, 3.10083e+08, 1.5747e+09, -1.09239e+09,
		-2.02373e+09, -1.09152e+09, -5.80176e+09, 0, -6.06332e+09, 0,
		-1.23945e+09, -9.46405e+08, -9.78925e+08, -1.00609e+09, -1.02956e+09,
		-1.53535e+09, -1.58479e+09, -1.63042e+09, -1.67294e+09, -1.43198e+09,
		-1.47061e+09, -1.5069e+09, -1.54113e+09, -1.5967e+09, -1.62614e+09,
		-1.56594e+09, -1.63982e+09, -1.76983e+09, -1.82778e+09, -1.58422e+09,
		-1.62898e+09, -1.10759e+09, -1.12237e+09, 1.0715e+08, 6.8374e+07,
		-7.83256e+10, -1.11839e+09, -1.13537e+09, -6.49394e+08, -6.82191e+09,
		-145075, -7.84922e+08, -2.96612e+09, 4.44918e+09, -1.07615e+09,
		-1.20791e+09, -1.26677e+09, -6.63725e+09, 1.87448e+07, 6.56029e+07,
		-7.34831e+09, 1.01175e+08, 1.30296e+08, 1.88709e+09, 3.37907e+09,
		-1.16851e+09, -1.5448e+09, -8.38734e+08, -4.45178e+09, 2.38419e-07,
		-4.64795e+09, 0, -9.49679e+08, -6.91726e+08, -7.22189e+08, -7.46534e+08,
		-7.67009e+08, -1.17632e+09, -1.21358e+09, -1.24796e+09, -1.28e+09,
		-1.10473e+09, -1.13243e+09, -1.15885e+09, -1.18398e+09, -1.22374e+09,
		-1.24551e+09, -1.26646e+09, -1.20344e+09, -1.32485e+09, -1.20494e+09,
		-1.23883e+09, -8.50787e+08, -8.61872e+08, 0, -3.34101e+07, 1.02561e+08,
		-1.04513e+10, -1.05684e+08, 3.52072e+09, -4.82333e+09, 1.00516e+07,
		9.18909e+06, 2.9511e+08, -5.04714e+07, -3.50544e+07, -2.78165e+07,
		1.34749e+08, 4.58223e+09, 1.31374e+08, 2.44558e+09, -3.36971e+09,
		-1.2137e+09, -344437, -179225, -122742, -93995.5, -76503, -64702.4,
		9.15549e+07, 1.53006e+08, -1.28875e+08, -5.32835e+07, 1.6812e+07,
		1.07465e+08, 8.90886e+06, -3.70765e+07, -3.32701e+07, -3.16828e+07,
		-3.02607e+07, -2.89788e+07, 4.82981e+07, -5.16505e+09, 3.95165e+06,
		5.06119e+09, 4.95694e+09, -5.5655e+09, 7.96447e+06, -4.1946e+07,
		7.47576e+07, -5.4172e+07, -3.06597e+07, -9.21333e+07, -4.67383e+07, 0,
		4.74651e+07, 71588.6, -90155.6, -68841.9, -38721.2, -46299.7, -161125,
		-127517, 4.33046e+07, 3.78118e+07, 1.96167e+07, -7.70638e+06,
		-8.97861e+07, -1.02138e+08, -9.5647e+07, -9.82185e+07, -4.40726e+07,
		-3.4978e+07, -2.67281e+07, -3.19669e+07, 172663, 110510, -5.58569e+09,
		-2.45799e+07, 5.54913e+09, 2.50284e+09, -8.33726e+07, -5.80616e+07,
		-4.14456e+07, -4.5922e+07, 2.24666e+09, 1.16258e+09, 4.04745e+08,
		-2.84264e+08, 470307, 373592, -21405.9, 114138, 135792, 126800, 104529,
		42645.3, -8.9818e+07, -7.86368e+07, -7.01449e+07, -6.34638e+07,
		-5.35975e+07, -4.98424e+07, -4.66369e+07, -4.38661e+07, -7.07606e+07,
		-8.23121e+07, -7.85765e+07, -8.20489e+07, 13630.3, -1.33173e+10,
		-1.3321e+08, -1.09418e+08, -6.02003e+07, -4.28928e+07, -3.38975e+07,
		-5.73595e+09, -2.96709e+09, -2.02373e+09, -1.5448e+09, 72537.6, 145075,
		217613, 1.00466e+06, 770498, 1.27724e+06, 440787, 150401, 0, -93265.1,
		-9.33118e+07, -8.16499e+07, -7.27954e+07, -6.58305e+07, -5.55487e+07,
		-5.16368e+07, -4.82981e+07, -4.54129e+07, -4.06712e+07, -3.8697e+07,
		-3.69301e+07, -3.53388e+07, -6.77048e+09, -8.78048e+07, 2.34409e+09,
		-3.12201e+09, 4.10576e+06, 4.48062e+07, 2.27591e+08, -3.78286e+07,
		4.82093e+07, 3.07791e+09, 4.62756e+07, 1.63308e+09, -2.25188e+09,
		-8.38734e+08, 41175.4, 51677.2, 66308.8, 89801.6, 135527, 268804,
		6.93243e+07, 1.73825e+08, -5.38463e+07, 0, 8.29104e+07, 1.15892e+08,
		3.68391e+07, 3.68391e+07, -3.18648e+09, 1.45031e+06, 3.31995e+09,
		3.1363e+09, -3.441e+09, 2.90981e+06, -4.27525e+07, 5.62425e+07,
		-3.13868e+07, -1.56934e+07, 0, 1.56934e+07, 59388.4, -67668, -26299.8,
		43427.6, 55617.8, -122419, -96768.4, 3.0189e+07, 2.65353e+07,
		1.29316e+07, -7.84382e+06, -4.80215e+07, -3.731e+07, -2.09956e+07,
		138336, 90181.5, -3.50235e+09, -1.88998e+07, 3.49618e+09, 1.53014e+09,
		-6.34947e+07, -4.48062e+07, -3.87505e+07, 1.37332e+09, 7.15443e+08,
		2.4789e+08, -1.819e+08, 367039, 289487, -8750.75, 91321.6, 106865,
		99293.9, 81532.2, 33927.1, -6.84522e+07, -6.01783e+07, -5.38463e+07,
		-4.88528e+07, -4.14552e+07, -3.86308e+07, -3.6215e+07, -3.4123e+07,
		-3.68099e+07, -2.92352e+07, -1.67666e+07, 11073.7, -8.23431e+09,
		-9.81908e+07, -8.12207e+07, -4.56121e+07, -3.28355e+07, -3.47348e+09,
		-1.8087e+09, -1.23945e+09, -9.49679e+08, 57133.4, 114267, 171400,
		749158, 567893, 946015, 328713, 112698, 0, -70374.3, -6.96907e+07,
		-6.13144e+07, -5.48492e+07, -4.97462e+07, -4.21895e+07, -3.93054e+07,
		-3.68391e+07, -3.47039e+07, -3.11856e+07, -2.9717e+07, -2.84004e+07,
		-2.72127e+07, -6.2031e+09, -7.55041e+07, -320676, -1.95035e+06, 0,
		-2.12873e+06, 0, -444223, -1.25828e+08, -3.63351e+07, 2.50773e+08,
		-3.85092e+07, 2.01152e+08, -3.01431e+07, -2.66871e+07, 1.95626e+08,
		-2.03301e+07, 1.52693e+08, -1.54142e+07, -2.1559e+07, 1.43806e+08,
		2.69309e+09, 1.22277e+08, 1.45667e+09, 1.08228e+08, 1.03124e+09,
		9.80725e+07, 8.10391e+08, 7.43795e+07, 8.28706e+07, -1.45484e+08,
		5.61345e+07, 1.35882e+08, 2.05679e+08, -1.06097e+08, 4.45569e+07,
		1.08992e+08, 1.67684e+08, 3.53389e+07, 4.0164e+07, 6.21804e+07,
		2.72127e+07, -316794, -333115, -344437, -354872, -364583, 7.49328e+07,
		6.87477e+07, 6.69333e+07, 6.11588e+07, 5.53962e+07, 5.01679e+07,
		-6.52205e+09, -7.89608e+07, -167247, -1.01148e+06, 0, -1.10197e+06, 0,
		-229787, 9.43045e+07, 8.7061e+07, 7.25758e+07, 5.95665e+07, 2.7332e+08,
		-4.63914e+07, 2.18761e+08, -3.62357e+07, 4.62218e+07, 2.07793e+08,
		-2.25585e+07, 1.61883e+08, -1.70788e+07, -2.21406e+07, 1.00916e+08,
		2.76937e+09, 9.06005e+07, 1.49655e+09, 8.29874e+07, 1.05392e+09,
		7.70413e+07, 8.2533e+08, -1.61877e+08, -4.09022e+07, 1.44852e+08,
		2.18546e+08, -1.13557e+08, -2.90063e+07, 1.13776e+08, 1.74929e+08,
		3.69301e+07, 4.13669e+07, 2.84004e+07, 566405, -173515, -179225,
		-184493, -189399, 2.67334e+07, 7.51737e+07, -2.30086e+07, 7.64519e+07,
		6.86326e+07, 6.31117e+07, 5.6163e+07, -6.72887e+09, 1.12059e+08,
		1.31926e+09, 9.26438e+08, -114716, -691135, 0, -752004, 0, -156728,
		9.55788e+07, 8.93997e+07, 7.35802e+07, 1.55807e+08, 3.00096e+08,
		-5.71037e+07, 2.39639e+08, -4.45025e+07, 1.19319e+08, 2.20651e+08,
		-2.58019e+07, 1.71598e+08, -1.94965e+07, -2.26125e+07, 8.03497e+07,
		2.84941e+09, 7.40857e+07, 1.53052e+09, 6.91443e+07, 1.07374e+09,
		1.27876e+09, -2.33555e+08, -4.61292e+07, 6.40675e+07, 2.33867e+08,
		-1.21944e+08, -3.14491e+07, 4.82945e+07, 1.83164e+08, 3.8697e+07,
		4.28322e+07, 2.9717e+07, 28817.5, 146344, -118915, -122742, -126275,
		-129567, 8.41992e+07, -2.44919e+07, -2.35072e+07, 7.82927e+07,
		6.38945e+07, -6.82072e+09, 1.97528e+08, 1.35649e+09, 9.39709e+08,
		-87954.2, -528331, -2.91038e-11, -574281, 0, -119638, 9.70141e+07,
		9.15605e+07, 7.44855e+07, 2.52233e+08, 3.33518e+08, -7.22383e+07,
		2.65638e+08, -5.6162e+07, 1.92532e+08, 2.35159e+08, -2.99354e+07,
		1.82543e+08, -2.25749e+07, -2.30042e+07, 6.78568e+07, 2.91726e+09,
		6.35758e+07, 1.5602e+09, 1.64489e+09, 1.58479e+09, 1.09144e+09,
		1.3051e+09, -2.0857e+08, -5.24944e+07, 2.27385e+07, 1.69246e+08,
		-1.316e+08, -3.41284e+07, 5.07965e+07, 1.25721e+08, 4.06712e+07,
		4.45699e+07, 3.11856e+07, 12109, 175376, 314032, -91113.9, -93995.5,
		-96656.3, -99136.5, -2.6044e+07, -2.4928e+07, -2.39199e+07,
		-6.89738e+09, 2.01249e+08, 1.39118e+09, 9.51877e+08, -429378, 0,
		-466325, 0, -97113.2, -1.74121e+08, 3.76922e+08, -9.4734e+07,
		2.99314e+08, -7.34558e+07, -1.42952e+08, 2.51911e+08, -3.5165e+07,
		1.9516e+08, -2.6466e+07, 4.65984e+07, 4.28928e+07, 3.28355e+07,
		-2.33412e+07, 5.9317e+07, 2.97659e+09, 2.37034e+09, 2.31413e+09,
		1.58674e+09, 1.68397e+09, 1.63042e+09, 1.10753e+09, 1.33125e+09,
		-6.06558e+07, 7.73716e+07, 1.85988e+08, 2.30944e+08, -3.71587e+07,
		5.37845e+07, 1.33131e+08, 2.03229e+08, 5671.81, 195410, 26452.3,
		68864.3, -74190.8, -76503, -78638.7, -80630.1, -2.77077e+07,
		-2.64401e+07, -2.53023e+07, -2.42746e+07, -6.96374e+09, -2.10895e+08,
		1.80061e+09, -2.85938e+09, 1.23534e+09, -1.88127e+09, -362688,
		-2.91038e-11, -393607, 0, -81943.9, 0, -8.13477e+07, 2.71602e+08,
		-4.18731e+07, 2.0996e+08, -3.1451e+07, -2.95247e+07, -2.36386e+07,
		4.50048e+09, 4.44741e+09, 3.02965e+09, 2.43282e+09, 2.38217e+09,
		1.61086e+09, 1.72149e+09, 1.67294e+09, 1.12232e+09, -3.08028e+09,
		1.49031e+08, 2.07356e+08, 3.02225e+08, -1.56582e+08, 5.73576e+07,
		1.41766e+08, 2.15528e+08, 4.54129e+07, 4.89562e+07, 3.47039e+07,
		2605.09, 210516, 15301.8, 80782.5, 223772, -64702.4, -66487, -68151.4,
		-2.80783e+07, -2.67892e+07, -2.56322e+07, -2.45874e+07, -9.89434e+09,
		-1.02608e+08, 4.14067e+09, -1.60614e+08, 264079, -297208, -196780,
		-74190.8, 72537.6, -8.47686e+07, -4.755e+07, 4.55379e+06, 1.06854e+09,
		2.48314e+09, 11321.8, 344437, 10978.6, 184493, 10651.8, 129567,
		-7.26567e+07, -6.38636e+07, -5.71709e+07, -5.18944e+07, -4.39548e+07,
		-4.09269e+07, -3.8339e+07, -3.60995e+07, -1.0182e+10, -1.05684e+08,
		2.15954e+09, -1.54669e+08, 402306, -309583, -260969, -93995.5, 9372.62,
		80372.8, -8.72402e+07, -4.88725e+07, -1.4619e+08, 4.30128e+09,
		2.20139e+08, 5.91324e+08, 10435.7, 354872, 10174.3, 189399,
		-7.47239e+07, -6.56409e+07, -5.87301e+07, -5.32835e+07, -4.51915e+07,
		-4.20626e+07, -3.93892e+07, -1.04401e+10, 2.49769e+08, 1.48125e+09,
		-1.49249e+08, 604602, -311867, -386788, -126275, 9009.97, 93112.5,
		151126, -8.9526e+07, -5.00326e+07, -1.34749e+08, 4.44747e+09,
		-1.4176e+08, 2.24047e+09, 1.92664e+09, 9710.44, 364583, -7.66364e+07,
		-6.72857e+07, -6.01734e+07, -5.45696e+07, -4.62786e+07, -4.31172e+07,
		-1.0673e+10, 1.54881e+08, 1.23493e+09, 7.51611e+07, -1.4439e+08, 922133,
		-705016, -754652, -189399, 8660.39, 112196, 176997, 78721.4, -157443,
		-146873, -5.11169e+07, -1.25384e+08, 4.58222e+09, -1.31374e+08,
		2.31421e+09, -1.37539e+08, 1.53546e+09, -7.84223e+07, -6.8822e+07,
		-6.15219e+07, -5.57716e+07, -4.72665e+07, -1.56629e+10, 1.70279e+09,
		-1.15428e+08, -7.03943e+07, -4.15697e+08, 0, 3.61185e+09, -1.23463e+08,
		-1.97516e+08, -2.01262e+07, 3.66985e+07, -4.17238e+07, -3.06597e+07,
		-2.47955e+08, 1.51461e+07, -1.4902e+09, -1.16131e+09, -3.53347e+07,
		-1.88097e+08, 1.10182e+08, 6.9956e+07, 3.34101e+07, 0, -4.4909e+08, 0,
		-9.33118e+07, -6.24591e+07, -3.29459e+08, 1.49012e-08, -3.42261e+08, 0,
		-6.96907e+07, -6.27866e+06, -2.10389e+08, -6.82595e+07, -2.04159e+08,
		-1.62338e+07, -2.00158e+08, -2.34971e+07, -1.95781e+08, 3.51808e+07,
		-1.94684e+08, -7.26567e+07, -7.47239e+07, -7.66364e+07, -7.84223e+07,
		-1.87665e+08, -3.86523e+07, -6.32614e+07, -6.39988e+07, -1.48223e+10,
		1.11481e+09, -1.00227e+08, -6.19191e+07, -2.91529e+08, 1.45765e+08,
		1.80205e+09, -1.09896e+08, -1.82557e+08, -6.41078e+07, 3.5903e+07,
		-3.13828e+07, 0, -2.33397e+08, 7.27064e+07, -1.14538e+08, 3.76668e+09,
		-1.1838e+09, -2.85644e+07, 3.70765e+07, 1.04934e+08, 6.68202e+07,
		3.19669e+07, -3.93184e+08, 0, -8.16499e+07, -2.89708e+08, 0,
		-3.00891e+08, 0, -6.13144e+07, -4.80237e+06, -1.91436e+08, -9.06091e+06,
		-1.86129e+08, -1.59466e+07, -1.80075e+08, -1.87636e+07, -1.78343e+08,
		-2.12642e+07, -1.7715e+08, -6.38636e+07, -6.56409e+07, -6.72857e+07,
		-6.8822e+07, -3.40082e+07, -5.56809e+07, -5.63148e+07, 3.48139e+07,
		-1.42115e+10, 1.17555e+09, -8.75679e+07, -5.54202e+07, -1.95262e+08,
		1.95262e+08, 1.18679e+09, -9.83864e+07, -1.71923e+08, -7.25616e+06,
		3.50544e+07, -2.36756e+07, 3.06597e+07, -2.294e+08, 9.87871e+07,
		-1.0702e+08, 3.90498e+09, -1.04233e+08, 1.88738e+09, 1.0023e+08,
		6.39339e+07, -3.50725e+08, 0, -7.27954e+07, -2.59422e+08, -1.49012e-08,
		-2.69231e+08, 0, -5.48492e+07, -3.80296e+06, -1.76513e+08, -1.02988e+07,
		-1.68721e+08, -6.0105e+07, -1.66412e+08, -1.53826e+07, 3.25727e+07,
		-1.753e+07, 3.11605e+07, -5.71709e+07, -5.87301e+07, -6.01734e+07,
		-6.15219e+07, -2.61065e+07, -2.83837e+07, -5.04572e+07, 3.40669e+07,
		3.31382e+07, -1.37951e+10, -1.29547e+08, 1.04221e+09, -5.02548e+07,
		-1.17902e+08, 1.76852e+08, -8.85236e+07, -1.64396e+08, -3.0934e+06,
		2.95247e+07, -1.47409e+07, 6.13193e+07, -2.32127e+08, 9.45935e+07,
		-1.00615e+08, 4.03072e+09, -9.89396e+07, 1.96292e+09, -9.51958e+07,
		1.24806e+09, -2.00134e+07, 3.32701e+07, 9.59008e+07, -3.17319e+08,
		1.49012e-08, -6.58305e+07, -2.35524e+08, 0, -2.44264e+08, 0,
		-4.97462e+07, -5.95554e+06, -1.60374e+08, -8.51681e+06, 3.34712e+07,
		-1.08101e+07, 3.20554e+07, -5.64848e+07, 3.07378e+07, -5.18944e+07,
		-5.32835e+07, -5.45696e+07, -5.57716e+07, -2.21356e+07, -2.41339e+07,
		3.2434e+07, 3.16376e+07, -1.31067e+10, -1.36789e+08, 1.1212e+09,
		-4.26074e+07, -1.11114e+08, -8.80424e+07, 4.44237e+09, -9.02719e+07,
		2.86811e+07, -7.86669e+06, 2.90599e+07, -1.82065e+07, -8.97877e+07,
		2.09357e+09, -8.81258e+07, 1.35031e+09, -8.51326e+07, -6.01222e+07,
		3.02607e+07, -2.49406e+08, -1.99525e+08, -2.67987e+08, -1.49012e-08,
		-5.55487e+07, -1.20054e+08, -1.20054e+08, -2.07276e+08, -1.49012e-08,
		-4.21895e+07, -2.17347e+06, 3.15911e+07, -4.24115e+06, 3.07055e+07,
		-6.13584e+06, 2.96946e+07, -9.45091e+06, 2.77077e+07, -1.09059e+07,
		2.67892e+07, -4.39548e+07, -4.51915e+07, -4.62786e+07, -4.72665e+07,
		-1.66253e+07, 2.96533e+07, -1.29549e+10, -1.38774e+08, 1.15599e+09,
		-3.96893e+07, -1.14608e+08, -7.40319e+07, 2.33481e+09, -8.58754e+07,
		2.81564e+07, -5.31483e+06, 2.84585e+07, -5.77759e+07, -8.60243e+07,
		4.53264e+09, -8.47647e+07, 1.5854e+09, -8.26673e+07, -1.30838e+07,
		2.89788e+07, -2.32132e+08, -2.78559e+08, -2.49212e+08, 0, -5.16368e+07,
		-7.462e+07, -1.1193e+08, -1.93154e+08, 0, -3.93054e+07, -1.86447e+06,
		2.97266e+07, -3.65703e+06, 2.90063e+07, -6.84144e+06, 2.72869e+07,
		-8.24848e+06, 2.64401e+07, -9.54858e+06, 2.56322e+07, -4.09269e+07,
		-4.20626e+07, -4.31172e+07, 4.15911e+07, -4.40244e+07, -1.26933e+10,
		-1.40271e+08, 1.37577e+09, -3.71943e+07, -1.17006e+08, 1.34842e+08,
		1.61475e+09, -8.16164e+07, 2.75159e+07, -3.18959e+06, 2.78165e+07,
		-1.15727e+07, -8.22465e+07, 4.61867e+09, -8.23396e+07, 2.37867e+09,
		-8.01269e+07, -2.1738e+08, -2.60856e+08, -2.33184e+08, 0, -4.82981e+07,
		-3.49927e+07, -6.99855e+07, -1.81075e+08, -1.49012e-08, -3.68391e+07,
		-1.61898e+06, 2.81076e+07, -4.65371e+06, 2.67954e+07, -6.0114e+06,
		2.6044e+07, -7.27045e+06, 2.53023e+07, -8.4401e+06, 2.45874e+07,
		-3.8339e+07, -3.93892e+07, -9.82749e+07, -8.09657e+07, -4.0365e+07,
		-8.25084e+07, -1.24169e+10, -1.41475e+08, 1.39512e+09, -3.50344e+07,
		-1.18752e+08, 1.31838e+08, -7.7629e+07, 2.66871e+07, -1.42051e+06,
		2.36386e+07, -7.5219e+06, -7.88599e+07, 4.70092e+09, -7.9132e+07,
		2.42068e+09, -7.86897e+07, 1.64298e+09, -3.70765e+07, -2.45546e+08,
		-2.45546e+08, -2.19331e+08, 0, -4.54129e+07, -1.70615e+08, 0,
		-3.47039e+07, -2.80931e+06, 2.6197e+07, -4.11291e+06, 2.55817e+07,
		-5.3292e+06, 2.4928e+07, -6.46325e+06, 2.42746e+07, -3.60995e+07,
		-1.01117e+08, 1.24532e+08, 2.44862e+08, -3.79841e+07, -3.88351e+07,
		-1.20152e+10, -1.43408e+08, 1.4333e+09, 1.28006e+08, 1.28006e+08,
		1.18464e+08, 1.18464e+08, 4.85558e+09, -7.30258e+07, 2.30042e+07,
		-4.28268e+06, -7.35364e+07, 2.49982e+09, -7.34617e+07, 1.69637e+09,
		-7.29457e+07, -3.32701e+07, -2.20363e+08, 0, -1.96557e+08, -1.57246e+08,
		-4.06712e+07, -9.20249e+07, -9.20249e+07, -3.11856e+07, -1.12219e+06,
		2.43073e+07, -2.23332e+06, 2.39636e+07, -3.28825e+06, 2.35072e+07,
		-2.77077e+07, -2.80783e+07, 1.13516e+08, 1.13516e+08, 2.24173e+08,
		-3.40669e+07, -3.48139e+07, -1.04527e+10, 1.44227e+08, 1.45194e+09,
		1.26803e+08, 1.26803e+08, 1.1059e+08, 1.1059e+08, 2.53722e+09,
		-7.10806e+07, 2.26125e+07, -2.96924e+06, -7.04896e+07, 4.9286e+09,
		-7.11284e+07, 1.72167e+09, 1.30771e+09, -3.16828e+07, -2.09868e+08,
		3.4978e+07, -1.87073e+08, -2.24487e+08, -3.8697e+07, -5.84704e+07,
		-8.77056e+07, -2.9717e+07, -1.00828e+06, 2.3299e+07, -2.01187e+06,
		2.30086e+07, -2.6044e+07, -2.64401e+07, -2.67892e+07, 1.1647e+08,
		1.1647e+08, 1.07971e+08, 1.07971e+08, -3.2434e+07, -3.31382e+07,
		-8.6256e+09, 1.44004e+08, 1.47022e+09, 1.25924e+08, 1.25924e+08,
		1.04786e+08, 1.04786e+08, 1.74614e+09, 1.74614e+09, 2.21406e+07,
		-1.823e+06, -6.81618e+07, 4.99909e+09, -6.88159e+07, 2.57335e+09,
		1.32621e+09, -3.02607e+07, -2.00461e+08, 6.68202e+07, -1.78583e+08,
		-2.14299e+08, -3.69301e+07, -2.79443e+07, -5.58886e+07, -2.84004e+07,
		-911499, 2.23875e+07, -2.44919e+07, -2.4928e+07, -2.53023e+07,
		-2.56322e+07, 1.14932e+08, 1.14932e+08, 1.08311e+08, 1.08311e+08,
		2.07543e+08, -3.09715e+07, -3.16376e+07, -5.9965e+09, 1.43954e+08,
		1.48811e+09, 1.25295e+08, 1.25295e+08, 1.00365e+08, 1.00365e+08,
		1.34416e+09, 2.1559e+07, -828524, -6.60159e+07, 5.06725e+09, 2.6083e+09,
		2.6083e+09, 1.76983e+09, 1.76983e+09, -2.89788e+07, -1.59835e+08,
		-7.45058e-09, -2.05122e+08, -2.05122e+08, -3.53388e+07, -2.72127e+07,
		-2.30086e+07, -2.35072e+07, -2.39199e+07, -2.42746e+07, -2.45874e+07,
		1.13743e+08, 1.13743e+08, 1.06516e+08, 1.06516e+08, 1.02276e+08,
		1.02276e+08, -2.96533e+07, -3.02854e+07, -1.0753e+10, 9.64915e+06,
		6.05993e+09, -6.74546e+09, 1.05058e+08, 6.60906e+07, 3.23905e+08,
		-6.21619e+07, -3.5903e+07, -2.84585e+07, 1.25384e+08, 4.70761e+09,
		-4.99433e+09, -1.6305e+09, 1.1634e+08, 344437, -10435.7, 1.72782e+08,
		-1.45922e+08, -6.01734e+07, 8.89414e+06, 1.14514e+08, 5.16368e+07,
		-3.45975e+07, -184493, -126275, -96656.3, -78638.7, -66487, -4.0365e+07,
		-3.79841e+07, -3.40669e+07, -3.2434e+07, -3.09715e+07, -2.96533e+07,
		-1.1034e+10, 1.25172e+08, 1.10135e+10, -1.10083e+10, 2.00216e+08,
		1.23068e+08, 3.60125e+08, -7.86858e+07, -3.66985e+07, -2.90599e+07,
		-9.7443e+09, -2.44558e+09, 1.11826e+08, 1.44871e+09, 179225, -10174.3,
		-1.6824e+08, -6.8822e+07, 8.62307e+06, 7.1156e+07, 5.55487e+07,
		8.28224e+06, -9710.44, 354872, -129567, -99136.5, -80630.1, -68151.4,
		-4.40244e+07, -8.25084e+07, -3.88351e+07, -3.48139e+07, -3.31382e+07,
		-3.16376e+07, -3.02854e+07, -6.92956e+09, 3.71778e+06, 3.9762e+09,
		-4.33185e+09, 7.99749e+07, 8.96124e+07, 2.49139e+08, -4.64862e+07,
		4.43418e+07, 3.12225e+09, -3.30096e+09, -1.10759e+09, 3.95702e+07,
		281564, -4480.31, 1.34414e+08, -1.15859e+08, 0, 4.88528e+07,
		1.24366e+08, 3.93054e+07, 3.93054e+07, 43629.5, 138657, 53800.8,
		92218.9, 215971, -7.07573e+09, 9.53683e+07, 7.1465e+09, -7.0257e+09,
		1.55921e+08, 1.34419e+08, 2.76136e+08, -5.86174e+07, -6.35844e+09,
		-1.63308e+09, 3.79796e+07, 9.18842e+08, 147469, -4300.6, -1.32451e+08,
		-5.63148e+07, 3.38904e+06, 9.77055e+07, 4.21895e+07, 4.21895e+07,
		-4120.88, 286044, 45387.9, 94326.6, 55508, 70266.1, -1.07784e+10,
		-1.64884e+08, 2.68532e+09, 1.85094e+07, -2.46476e+08, -4.28314e+08,
		-6.15217e+08, -1.02862e+09, -7.5864e+07, 3.68522e+07, 4.23493e+07,
		3.94468e+07, -7.85897e+06, -1.30971e+08, -2.9562e+09, -2.81069e+09,
		-2.44334e+09, -1.86594e+09, -397204, -305253, -204906, -100152, 7389.5,
		-40404.5, 203614, -102800, 16388.8, -8.53144e+09, -2.79447e+07,
		-2.55312e+09, -2.57133e+09, 4.51397e+09, 2.87909e+07, -1.08207e+06,
		3.08343e+07, 4.67283e+07, 4.64663e+07, -2.7486e+06, -5.96413e+06,
		7.51645e+07, -5.22096e+09, -2.39763e+09, -1.29654e+09, -3.20344e+08,
		3.23801e+07, -8.44695e+06, -3.92225e+07, -1.6659e+07, -251752, -285572,
		-255681, -195734, -112688, -78721.4, 80437, -3.37339e+06, -4.1427e+06,
		-5.2155e+06, -6.76612e+06, -1.13262e+10, -1.12248e+08, -2.71487e+09,
		2.78275e+09, -2.80714e+09, 7.96447e+06, -6.60261e+07, 8.20425e+06,
		-4.64331e+07, -3.06597e+07, -5.12883e+09, -3.94152e+09, -2.64581e+09,
		-1.2932e+09, -27477.7, -189226, -126285, -61885.7, -52769.9, -228944,
		-217675, -4.65134e+07, -9.6172e+07, -1.14189e+08, -1.03612e+08,
		-6.98336e+07, -7.42824e+07, -6.95615e+07, -8.18487e+07, -3.67272e+07,
		-3.4978e+07, -4.00921e+07, -3.19669e+07, -65816.8, -28364.8,
		-1.2301e+10, -2.18209e+08, -2.94107e+09, 1.25421e+09, -1.2714e+08,
		-5.80616e+07, -4.14456e+07, -3.93617e+07, -3.2507e+09, -3.6874e+09,
		-3.30144e+09, -2.52738e+09, -186796, -93611.3, -404340, -185685,
		-100411, -24809.2, 67222.6, 29892.3, -8.9818e+07, -7.86368e+07,
		-7.01449e+07, -6.34638e+07, -5.35975e+07, -4.98424e+07, -4.66369e+07,
		-4.38661e+07, -5.5036e+07, -5.98633e+07, -5.71465e+07, -6.8374e+07,
		-13630.3, -7.18057e+09, -9.05436e+07, -1.73766e+09, 1.7205e+09,
		-1.71553e+09, 2.90981e+06, -4.37057e+07, 5.97253e+06, -3.26758e+09,
		-2.54293e+09, -1.73007e+09, -8.66817e+08, -22689.3, -149202, -124650,
		-108645, -77857.9, -178334, -170770, -3.82631e+07, -7.58814e+07,
		-8.95769e+07, -8.1391e+07, -3.20144e+07, -2.9848e+07, -2.09956e+07,
		-52132.3, -22847.3, -7.6826e+09, -1.63971e+08, -1.80792e+09,
		7.54534e+08, -9.5983e+07, -4.48062e+07, -2.58337e+07, -2.0198e+09,
		-2.29506e+09, -2.06153e+09, -1.58231e+09, -144743, -70502.1, -306587,
		-142816, -76944.2, -18007.4, 53382.5, 23206.3, -6.84522e+07,
		-6.01783e+07, -5.38463e+07, -4.88528e+07, -4.14552e+07, -3.86308e+07,
		-3.6215e+07, -3.4123e+07, -2.454e+07, -2.33882e+07, -1.67666e+07,
		-11073.7};
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
				XOLOTL_REQUIRE_CLOSE_ZT(knownPartials[startingIdx + j],
					hPartials[startingIdx + index], 0.01, 1.0e-6);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 21);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.14068, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 1.9, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 29106.01, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 103680000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto& reg0 = cluster.getRegion();
	Composition lo = reg0.getOrigin();
	Composition hi = reg0.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::V] = 0;
	comp[Spec::I] = 2;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 10);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.17725, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 5.4559, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 156594191.2, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 51840000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.5, 0.01);
	const auto& reg1 = cluster.getRegion();
	lo = reg1.getOrigin();
	hi = reg1.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::I] = 0;
	comp[Spec::Void] = 15;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 19);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.361465, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 20.96172, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg2 = cluster.getRegion();
	lo = reg2.getOrigin();
	hi = reg2.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 15);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 20);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Void] = 0;
	comp[Spec::Frank] = 12;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.60856, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 11.86668, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg3 = cluster.getRegion();
	lo = reg3.getOrigin();
	hi = reg3.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 10);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 15);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Frank] = 0;
	comp[Spec::Perfect] = 5;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 5);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.321142, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.2480, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 6150067.77, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 20736000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.7, 0.01);
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 6);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::Perfect] = 0;
	comp[Spec::Faulted] = 20;
	cluster = network.findCluster(comp, plsm::onHost);
	BOOST_REQUIRE_EQUAL(cluster.getId(), 16);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.78702848, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 14.63612, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg5 = cluster.getRegion();
	lo = reg5.getOrigin();
	hi = reg5.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Void], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Void], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Faulted], 20);
	BOOST_REQUIRE_EQUAL(hi[Spec::Faulted], 21);
	BOOST_REQUIRE_EQUAL(lo[Spec::Frank], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Frank], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::Perfect], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::Perfect], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
}

BOOST_AUTO_TEST_SUITE_END()
