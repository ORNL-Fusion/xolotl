#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <xolotl/core/network/AlloyReactionNetwork.h>

using namespace std;
using namespace xolotl::core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the Alloy network.
 */
BOOST_AUTO_TEST_SUITE(AlloyNetwork_testSuite)

BOOST_AUTO_TEST_CASE(fullyRefined) {

	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 6 4" << std::endl << "process=reaction"
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

	using NetworkType = AlloyReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxSize = opts.getMaxImpurity();
	NetworkType network( { maxV, maxSize, maxSize, maxI, maxSize, maxSize }, 1,
			opts);

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 23);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 23);
	// TODO: check it is within a given range?
	auto deviceMemorySize = network.getDeviceMemorySize();
	BOOST_REQUIRE(deviceMemorySize > 140000);
	BOOST_REQUIRE(deviceMemorySize < 160000);

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
	knownDFill[0] = { 0, 9, 10, 11, 17, 18, 19, 20, 1, 2, 3, 12 };
	knownDFill[1] = { 1, 0, 9, 10, 17, 18, 19, 20, 21, 2, 3, 12, 11 };
	knownDFill[2] = { 2, 0, 10, 1, 9, 17, 18, 19, 20, 21, 22, 3, 11, 12 };
	knownDFill[3] = { 3, 0, 11, 1, 10, 2, 9, 17, 18, 19, 20, 21, 22, 4, 12 };
	knownDFill[4] = { 4, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 5, 6, 7, 8,
			20 };
	knownDFill[5] =
			{ 5, 4, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 6, 7, 8 };
	knownDFill[6] =
			{ 6, 4, 10, 5, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 7, 8 };
	knownDFill[7] = { 7, 4, 11, 5, 10, 6, 9, 13, 14, 16, 17, 18, 19, 20, 21, 22,
			8 };
	knownDFill[8] = { 8, 4, 12, 5, 11, 6, 10, 7, 9, 13, 14, 17, 18, 19, 20, 21,
			22 };
	knownDFill[9] = { 9, 0, 20, 1, 21, 2, 22, 4, 19, 5, 6, 7, 8, 13, 10, 11, 12,
			14, 15, 16, 17, 18 };
	knownDFill[10] = { 10, 0, 19, 1, 20, 2, 21, 3, 22, 4, 18, 5, 6, 7, 8, 9, 11,
			12, 13, 14, 15, 16, 17 };
	knownDFill[11] = { 11, 0, 18, 1, 19, 2, 20, 3, 21, 4, 17, 5, 6, 7, 8, 9, 10,
			12, 13, 14, 15, 16, 22 };
	knownDFill[12] = { 12, 0, 17, 1, 18, 2, 19, 3, 20, 4, 5, 6, 7, 8, 9, 11, 10,
			13, 14, 15, 16, 21, 22 };
	knownDFill[13] = { 13, 4, 5, 6, 7, 8, 9, 14, 10, 11, 12, 17 };
	knownDFill[14] = { 14, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17 };
	knownDFill[15] = { 15, 4, 5, 6, 9, 16, 10, 11, 12, 17, 22, 18, 21, 19, 20 };
	knownDFill[16] =
			{ 16, 4, 5, 6, 7, 9, 10, 11, 12, 15, 17, 18, 22, 19, 21, 20 };
	knownDFill[17] = { 17, 0, 1, 2, 3, 4, 5, 6, 13, 15, 7, 14, 16, 8, 9, 18, 10,
			19, 11, 20, 12, 21, 22 };
	knownDFill[18] = { 18, 0, 1, 2, 3, 4, 5, 13, 15, 6, 14, 16, 7, 8, 9, 19, 10,
			20, 11, 21, 12, 22, 17 };
	knownDFill[19] = { 19, 0, 1, 2, 3, 4, 13, 15, 5, 14, 16, 6, 7, 8, 9, 20, 10,
			21, 11, 22, 12, 17, 18 };
	knownDFill[20] = { 20, 0, 1, 2, 3, 4, 14, 16, 5, 6, 7, 8, 9, 21, 10, 22, 11,
			13, 15, 12, 17, 19, 18 };
	knownDFill[21] = { 21, 1, 2, 3, 6, 7, 8, 9, 22, 10, 13, 15, 11, 14, 16, 12,
			17, 20, 18, 19 };
	knownDFill[22] = { 22, 2, 3, 7, 8, 9, 13, 15, 10, 14, 16, 11, 12, 17, 21,
			18, 20, 19 };
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 413);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto &row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = { 1000.0 };
	network.setTemperatures(temperatures);
	network.syncClusterDataOnHost();
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 5048430879.96, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
	Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(network.getTotalConcentration(dConcs, Spec::Void, 1),
			2.0, 0.01);
	BOOST_REQUIRE_CLOSE(
			network.getTotalAtomConcentration(dConcs, Spec::Perfect, 3), 30.0,
			0.01);
	BOOST_REQUIRE_CLOSE(
			network.getTotalRadiusConcentration(dConcs, Spec::Frank, 6),
			1.394446, 0.01);
	BOOST_REQUIRE_CLOSE(
			network.getTotalTrappedAtomConcentration(dConcs, Spec::Faulted, 1),
			0.0, 0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(network.getLeftSideRate(dConcs, 0, gridId), 7782153079,
			0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
	Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = { -1.79219e+09, 1.85298e+09, 3.81785e+09,
			9.93479e+09, -9.45535e+09, -4.17172e+09, -7.68927e+08, 2.88205e+09,
			8.58606e+09, -7.85155e+10, -3.45865e+10, -2.27619e+10, -1.39981e+10,
			-5.3704e+09, -1.02448e+10, -3.69144e+09, -6.92805e+09, 4.06728e+08,
			2.84616e+08, 1.50115e+09, 3.00266e+09, 3.1405e+09, 4.93097e+09, 0 };
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials =
			{ -7.78215e+09, -5.78851e+08, 2.67961e+08, 9.4625e+08, 12438.8,
					11978.4, 11558.6, -87954.2, 333115, 179225, 126275,
					3.56182e+09, -6.54255e+09, 4.14067e+09, -1.60614e+08,
					4.55379e+06, 11321.8, 10978.6, -118915, -91113.9, -74190.8,
					344437, 184493, 2.24503e+09, 4.01861e+09, -4.44836e+09,
					2.15954e+09, 2.15954e+09, 4.30128e+09, -1.4619e+08, 10435.7,
					-179225, -122742, -93995.5, -76503, -64702.4, 354872,
					1.80502e+09, 1.80502e+09, -907422, 1.48125e+09, 1.48125e+09,
					2.24047e+09, 2.24047e+09, 4.44747e+09, 4.44747e+09, -354872,
					-184493, -126275, -96656.3, -78638.7, -66487, 3.58296e+08,
					3.17469e+09, -1.03899e+10, -3.76754e+09, -2.02102e+09,
					-1.42287e+09, -1.11839e+09, -1.05684e+08, -1.08527e+08,
					-8.78048e+07, -8.90013e+07, -1.25903e+07, -2.23786e+07,
					-3.00636e+07, 6.29138e+07, 5.65822e+07, 5.13656e+07,
					4.71167e+07, 4.71167e+07, -8.08967e+09, 3.76754e+09,
					-1.34104e+08, -2.08005e+09, -1.45647e+09, -8.72402e+07,
					-8.9526e+07, -7.27232e+07, -7.36831e+07, -7.15151e+07,
					-1.56958e+07, -2.1501e+07, -6.9211e+07, 5.43125e+07,
					4.99631e+07, 4.61292e+07, -6.78453e+09, 2.02102e+09,
					-1.15428e+08, 3.90164e+09, -1.23463e+08, -7.47239e+07,
					-7.66364e+07, -6.24591e+07, -6.32614e+07, -6.27866e+06,
					-6.82595e+07, -5.82218e+07, -5.9539e+07, -6.06558e+07,
					4.80338e+07, 4.49048e+07, -4.68922e+09, 1.42287e+09,
					1.42287e+09, 2.08005e+09, 2.08005e+09, 4.02511e+09,
					-1.14538e+08, -6.56409e+07, -6.72857e+07, -5.56809e+07,
					-4.80237e+06, -4.99631e+07, -1.02731e+08, -5.24944e+07,
					-5.34523e+07, -5.42917e+07, 4.32314e+07, -4.44048e+08,
					1.47669e+09, 1.11839e+09, 1.45647e+09, 1.45647e+09,
					2.13645e+09, 2.13645e+09, 4.13965e+09, 4.13965e+09,
					-5.87301e+07, -6.01734e+07, -4.32314e+07, -4.49048e+07,
					-4.61292e+07, -9.42333e+07, -4.79553e+07, -4.86906e+07,
					-9.2016e+10, -4.14058e+09, -2.78022e+09, -4.30121e+09,
					-2.85665e+09, -4.44741e+09, -2.92234e+09, -3.68611e+09,
					-1.91407e+09, -3.83243e+09, -3.96445e+09, -4.08535e+09,
					5.87301e+07, -4.38874e+09, -2.66788e+09, -2.72098e+09,
					-2.78814e+09, -4.58222e+09, -3.0297e+09, -3.07791e+09,
					-1.19014e+09, -1.69712e+09, -4.72394e+10, -2.15943e+09,
					-1.38885e+09, -2.24038e+09, -1.43713e+09, 76503,
					-1.47715e+09, 66487, -1.51157e+09, -1.94206e+09,
					-5.78712e+08, -2.01242e+09, -2.07691e+09, 5.34523e+07,
					4.86906e+07, 6.09466e+09, -1.50422e+09, -1.49655e+09,
					-2.31421e+09, -2.38223e+09, -1.58681e+09, -1.61091e+09,
					-4.11028e+08, -3.17956e+10, -1.48108e+09, -9.0531e+08,
					118915, -9.46405e+08, 93995.5, -9.78925e+08, 78638.7,
					-1.00609e+09, -1.34737e+09, -1.35153e+08, -1.39081e+09,
					5.82218e+07, 5.24944e+07, 4.79553e+07, 3.10083e+08,
					1.5747e+09, -1.09239e+09, -1.58488e+09, -1.6305e+09,
					-1.09152e+09, -1.10759e+09, -1.07389e+09, -2.38957e+10,
					320676, -6.49394e+08, 173515, -6.91726e+08, 122742,
					-7.22189e+08, 96656.3, -7.46534e+08, -1.11839e+09,
					6.29138e+07, 5.65822e+07, 5.13656e+07, 4.71167e+07,
					1.30296e+08, 1.88709e+09, 3.37907e+09, -1.2137e+09,
					-1.24806e+09, -8.38734e+08, -8.50787e+08, -8.10701e+08,
					-8.25523e+08, -9.95263e+09, -1.05684e+08, -8.72402e+07,
					-7.47239e+07, -6.56409e+07, -5.87301e+07, 1.34749e+08,
					4.58223e+09, -2.31421e+09, -1.58488e+09, -1.2137e+09,
					-344437, -1.02452e+10, -1.08527e+08, -8.9526e+07,
					-7.66364e+07, -6.72857e+07, -6.01734e+07, -4.58222e+09,
					-2.38223e+09, -1.6305e+09, -1.24806e+09, 344437, 344437,
					-6.77003e+09, -8.78048e+07, -7.27232e+07, -6.24591e+07,
					4.82093e+07, 3.07791e+09, -1.58681e+09, -1.09152e+09,
					-8.38734e+08, 41175.4, 322739, 199146, 199146, 167749,
					167749, -6.92883e+09, -8.90013e+07, -7.36831e+07,
					-6.32614e+07, -5.56809e+07, -3.07791e+09, -1.61091e+09,
					-1.10759e+09, -8.50787e+08, 281564, 281564, 193339, 193339,
					156734, 156734, 295039, -5.84714e+09, -320676, -333115,
					-344437, -354872, -7.55041e+07, -1.25828e+08, 8.28706e+07,
					7.43795e+07, 6.21804e+07, 7.49328e+07, 6.72882e+07,
					5.56823e+07, -4.32314e+07, 1.43806e+08, 2.69309e+09,
					1.22277e+08, 1.45667e+09, 1.08228e+08, 1.03124e+09,
					9.80725e+07, 8.10391e+08, -316794, -6.16223e+09, -167247,
					-173515, -179225, -184493, -7.89608e+07, 9.43045e+07,
					8.72402e+07, 7.27232e+07, 2.67334e+07, 7.66364e+07,
					6.32614e+07, -4.99631e+07, -4.49048e+07, 1.00916e+08,
					2.76937e+09, 9.06005e+07, 1.49655e+09, 8.29874e+07,
					1.05392e+09, 7.70413e+07, 8.2533e+08, 566405, -6.36287e+09,
					-114716, -118915, -122742, -126275, 1.12059e+08,
					1.31938e+09, 9.26539e+08, 9.55788e+07, 8.9526e+07,
					7.36831e+07, -5.82218e+07, -1.02731e+08, -4.61292e+07,
					8.03497e+07, 2.84941e+09, 7.40857e+07, 1.53052e+09,
					6.91443e+07, 1.07389e+09, 1.27876e+09, 28817.5, 146344,
					-6.44874e+09, -87954.2, -91113.9, -93995.5, -96656.3,
					1.97528e+08, 1.35658e+09, 9.39788e+08, -6.9211e+07,
					-5.9539e+07, -5.24944e+07, -9.42333e+07, 6.78568e+07,
					2.91739e+09, 6.35758e+07, 1.56033e+09, 1.64489e+09,
					1.58488e+09, 1.09152e+09, 1.3051e+09, 12109, 175376, 314032,
					-6.47587e+09, -74190.8, -76503, -78638.7, -6.06558e+07,
					-5.34523e+07, -4.79553e+07, 5.9317e+07, 2.9767e+09,
					2.37034e+09, 2.31421e+09, 1.58681e+09, 1.68397e+09,
					1.6305e+09, 1.10759e+09, -8.10701e+08, 5671.81, 331016,
					26452.3, 68864.3, -6.54007e+09, -64702.4, -66487,
					-5.42917e+07, -4.86906e+07, 4.50048e+09, 4.44747e+09,
					3.0297e+09, 2.43282e+09, 2.38223e+09, 1.61091e+09,
					-1.07389e+09, -8.25523e+08, 2605.09, 325344, 15301.8,
					208641, 373467 };
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto &row = rowIter->second;
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
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 29106.01,
			0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 103680000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto &reg0 = cluster.getRegion();
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
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 156594191.2,
			0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 51840000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.5, 0.01);
	const auto &reg1 = cluster.getRegion();
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
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0,
			0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(cluster.getMigrationEnergy(),
			std::numeric_limits<double>::infinity());
	const auto &reg2 = cluster.getRegion();
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
	BOOST_REQUIRE_EQUAL(cluster.getMigrationEnergy(),
			std::numeric_limits<double>::infinity());
	const auto &reg3 = cluster.getRegion();
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
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 6150067.77, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 20736000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.7, 0.01);
	const auto &reg4 = cluster.getRegion();
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
	BOOST_REQUIRE_EQUAL(cluster.getMigrationEnergy(),
			std::numeric_limits<double>::infinity());
	const auto &reg5 = cluster.getRegion();
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

BOOST_AUTO_TEST_SUITE_END()
