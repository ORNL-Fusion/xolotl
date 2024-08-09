#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/options/ConfOptions.h>
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
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=5 0 0 5 1" << std::endl
			  << "reactionFilePath=reactionRates.txt" << std::endl
			  << "process=reaction" << std::endl
			  << "tempParam=1500" << std::endl;
	paramFile.close();

	// Create the reaction file
	std::string reactionFile = "reactionRates.txt";
	std::ofstream reacFile(reactionFile);
	reacFile << "0 0 1 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
			 << std::endl
			 << "0 1 0 -0.160789024638881 -12.7387253625273 4.2488E+11 4.23317"
			 << std::endl
			 << "0 2 0 2.87060634573234 -8.37649260139129 9.9337E+11 3.43830674"
			 << std::endl
			 << "1 0 0 9.18245700000039 3.7599999999763 2.9921E+11 1"
			 << std::endl
			 << "1 1 0 3.14296517286615 -1.4324963007117 0 0" << std::endl
			 << "1 2 0 2.26277914822638 -9.53632166324027 6.4829E+11 3.69567683"
			 << std::endl
			 << "1 3 0 4.17605871610343 -5.00715984031393 7.4802E+10 5.34571352"
			 << std::endl
			 << "Reactions" << std::endl
			 << "0 0 1 0.0001" << std::endl
			 << "0 1 0 0.0001" << std::endl
			 << "0 2 0 0.0001" << std::endl;
	reacFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create the network
	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxXe, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{maxXe + 1, maxV + 1, maxI + 1}};
	NetworkType network(maxSpeciesAmounts, subdivRatios, 1, opts);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 29);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 29);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.547, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0409168, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 29);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 2);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 20,
		21, 22, 23, 25, 26, 27, 28, 3, 9, 14, 19, 24};
	knownDFill[1] = {1, 0, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 17, 19,
		20, 21, 22, 24, 25, 26, 27, 8, 13, 18, 23, 28};
	knownDFill[2] = {2, 0, 1, 3, 4, 5, 6, 9, 10, 11, 14, 15, 16, 19, 20, 21, 24,
		25, 26, 7, 8, 12, 13, 17, 18, 22, 23, 27, 28};
	knownDFill[3] = {3, 0, 4, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	knownDFill[4] = {4, 0, 5, 1, 3, 2, 6, 9, 11, 12};
	knownDFill[5] = {5, 0, 6, 1, 4, 2, 3, 9, 10, 11, 14, 15, 16, 19, 20, 21, 7};
	knownDFill[6] = {6, 0, 7, 1, 5, 2, 4, 3, 9, 10, 14, 15, 19, 20, 8, 11, 12};
	knownDFill[7] = {7, 0, 8, 1, 6, 2, 5, 3, 12};
	knownDFill[8] = {8, 0, 1, 7, 2, 6, 3, 13};
	knownDFill[9] = {9, 0, 10, 1, 2, 3, 4, 5, 6, 11, 14};
	knownDFill[10] = {10, 0, 11, 1, 9, 2, 3, 5, 6, 12, 15};
	knownDFill[11] = {11, 0, 12, 1, 10, 2, 9, 3, 6, 4, 5, 13, 16};
	knownDFill[12] = {12, 0, 13, 1, 11, 2, 10, 3, 7, 4, 6, 5, 17};
	knownDFill[13] = {13, 0, 1, 12, 2, 11, 3, 8, 5, 6, 18};
	knownDFill[14] = {14, 0, 15, 1, 2, 3, 9, 5, 6, 16, 19};
	knownDFill[15] = {15, 0, 16, 1, 14, 2, 3, 10, 5, 6, 17, 20};
	knownDFill[16] = {16, 0, 17, 1, 15, 2, 14, 3, 11, 5, 9, 18, 21};
	knownDFill[17] = {17, 0, 18, 1, 16, 2, 15, 3, 12, 5, 10, 6, 9, 22};
	knownDFill[18] = {18, 0, 1, 17, 2, 16, 3, 13, 5, 11, 6, 10, 23};
	knownDFill[19] = {19, 0, 20, 1, 2, 3, 14, 5, 6, 21, 24};
	knownDFill[20] = {20, 0, 21, 1, 19, 2, 3, 15, 5, 6, 22, 25};
	knownDFill[21] = {21, 0, 22, 1, 20, 2, 19, 3, 16, 5, 14, 23, 26};
	knownDFill[22] = {22, 0, 23, 1, 21, 2, 20, 3, 17, 5, 15, 6, 14, 27};
	knownDFill[23] = {23, 0, 1, 22, 2, 21, 3, 18, 5, 16, 6, 15, 28};
	knownDFill[24] = {24, 0, 25, 1, 2, 3, 19, 26};
	knownDFill[25] = {25, 0, 26, 1, 24, 2, 3, 20, 27};
	knownDFill[26] = {26, 0, 27, 1, 25, 2, 24, 3, 21, 5, 19, 28};
	knownDFill[27] = {27, 0, 28, 1, 26, 2, 25, 3, 22, 5, 20, 6, 19};
	knownDFill[28] = {28, 0, 1, 27, 2, 26, 3, 23, 5, 21, 6, 20};
	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 416);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1500.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 6821697397703317, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::Xe, 1), 26.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::Xe, 3), 60.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1), 12.106857,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Xe, 1), 75.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 2704.45, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {2.36776e+06, -1.1276e+11, -1.1276e+11,
		2.57526e+16, 1.44634e+11, -9.12248e+09, -1.28282e+11, -1.33376e+11,
		-1.37678e+11, 6.31097e+15, 385105, -360.232, 3454.38, 3321.87, -390996,
		1.56336e+11, -360.231, 3178.43, 3321.87, -390996, 6.82154e+15, 528872,
		3178.42, 3321.87, -6.31101e+15, -6.82158e+15, 1.28282e+11, 1.33376e+11,
		1.37678e+11, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-2704.45, -8.50534, -2317.66, -6.39485,
		-230.224, -7.31052, -7.61095, -7.86472, -6.93297, -7.31045, -7.61095,
		-7.86472, -6.93296, -7.31045, -7.61095, -7.86472, 409154, -7.31045,
		-7.61095, -7.86472, 409154, -2.91825, -7.61095, -7.86472, 2.48123e-15,
		388034, 388034, 388034, 388034, -1.1276e+11, 2309.16, 2317.67,
		-1.1276e+11, -2.47254, -261.278, -2.70373, -2.87317, -2.47259, -2.65018,
		-2.77463, -2.87338, -2.47259, -2.65018, -2.77476, -2.8739, -2.47259,
		-2.65018, -2.77476, -2.87394, -2.47259, -2.65018, -2.77476, -2.87394,
		0.000781093, 0.000738891, 0.000305525, 1.13101e-05, 7.41541e-09,
		-1.1276e+11, -2317.66, 9.89035, -1.1276e+11, -2707.82, -3160.94,
		-3038.74, -2707.82, -2902.31, -3038.74, -2707.82, -2902.31, -3038.74,
		-2707.82, -2902.31, -3038.74, -2707.82, -2902.31, -3038.74, 4.9183e-09,
		4.61987e-09, 2.36938e-09, 4.13414e-09, 9.24533e-12, 7.06832e-10,
		1.12611e-12, 1.17904e-12, 1.12611e-12, 1.17904e-12, -2.76144e+12,
		6.39485, -1.1276e+11, -1.1276e+11, -1.1276e+11, -1.21883e+11,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 3.18738e+10, -1.21883e+11,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 6.31101e+15, -1.21883e+11,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 6.31101e+15, 3.44525e+10,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 6.31112e+15, 6.8217e+15,
		120071, 5.17505e-05, 5.34197e-05, -1.1276e+11, 223.83, -28.4037,
		1.1276e+11, 2.19519, -2707.82, -8.93039e-05, 1.44634e+11, 1.00349e-13,
		3.46518e-20, -1.21883e+11, -222.914, -280.736, -258.806, -256.156,
		1.1276e+11, -9.12247e+09, -258.629, -275.96, -288.117, -258.629,
		-275.96, -288.117, -258.629, -275.96, -288.117, 4.9183e-09,
		-1.28282e+11, 0.300431, 7.61172, 258.504, -26.8386, -330.92, 2707.82,
		-1.28282e+11, -8.9308e-05, -9.5024e-05, -8.9308e-05, -9.5024e-05,
		-8.9308e-05, -9.5024e-05, 4.61987e-09, 4.97738e-05, 3.46518e-20,
		-1.33376e+11, 0.253768, 7.8655, -0.0990872, 2.77485, 3160.94, 3160.94,
		-1.33376e+11, 5.17505e-05, -1.37678e+11, -7.86472, 2.87394, 2.87394,
		3038.74, 3038.74, -1.37678e+11, 5.34197e-05, -2.57395e+11, 6.93297,
		6.93297, -2.47259, -2707.82, 0, 1.1276e+11, -258.629, -8.9308e-05,
		1.43843e-10, 6.31112e+15, -1.21883e+11, 0.377476, 7.31058, -0.177594,
		388036, -2902.31, 231.557, 1.21883e+11, -9.5024e-05, 2.36938e-09,
		17.7205, -1.28282e+11, 0.300508, 7.61151, -0.124578, 2.65018, -330.92,
		2707.82, 7.62939e-05, 1.28282e+11, 258.629, -29.4887, 4.13414e-09,
		4.97738e-05, -1.33376e+11, 0.253768, 7.86546, -0.0991765, 2.77476,
		2902.31, 2902.31, 0, 1.33376e+11, 8.9308e-05, 8.9308e-05, 1103.84,
		5.17505e-05, -1.37678e+11, -7.86472, 2.87394, 2.87394, 3038.74, 3038.74,
		0, 1.37678e+11, 288.117, 288.117, 5.34197e-05, -6.31123e+15, 6.93297,
		6.93297, -2.47259, -2707.82, 0, 1.1276e+11, -258.629, -8.9308e-05,
		1.05064e-12, 6.31112e+15, -1.21883e+11, 0.377476, 7.31045, -0.177594,
		388036, -2902.31, 0, 1.21883e+11, -275.96, -9.5024e-05, 9.24533e-12,
		1.56335e+11, -1.28282e+11, 0.300508, 7.61099, -0.124578, 2.65019,
		-330.92, 2707.82, 0, 1.28282e+11, -29.4887, 258.629, 7.06832e-10,
		0.00120217, -1.33376e+11, 0.253768, 7.86503, -0.0991765, 2.77476,
		2902.31, 2902.31, 0, 1.33376e+11, 275.96, 275.96, 8.9308e-05,
		8.9308e-05, 5.17505e-05, -1.37678e+11, -7.86472, 2.87394, 2.87394,
		3038.74, 3038.74, 0, 1.37678e+11, 288.117, 288.117, 9.5024e-05,
		9.5024e-05, 5.34197e-05, -6.31123e+15, 6.93297, 6.93297, -2.47259,
		-2707.82, 0, 1.1276e+11, -258.629, -8.9308e-05, 1.05064e-12,
		6.31112e+15, -2.78218e+11, 0.377476, 7.31045, -0.177594, 388036,
		-2902.31, 0, 1.21883e+11, -275.96, -9.5024e-05, 1.12611e-12, 6.8217e+15,
		-1.28282e+11, 0.300508, 7.61095, -0.124578, 409164, -330.92, 2707.82, 0,
		1.28282e+11, -29.4887, 258.629, 1.17904e-12, 120071, -1.33376e+11,
		0.253768, 7.86473, -0.0991765, 2.77476, 2902.31, 2902.31, 0,
		1.33376e+11, 275.96, 275.96, 8.9308e-05, 8.9308e-05, 5.17505e-05,
		-1.37678e+11, -7.86472, 2.87394, 2.87394, 3038.74, 3038.74, 0,
		1.37678e+11, 288.117, 288.117, 9.5024e-05, 9.5024e-05, 5.34197e-05,
		-6.31112e+15, 6.93297, 6.93297, -2.47259, -2707.82, 1.1276e+11,
		1.1276e+11, 1.05064e-12, -6.8217e+15, 0.377476, 7.31045, -0.177594,
		388036, -2902.31, 1.21883e+11, 1.21883e+11, 1.12611e-12, -123124,
		0.300508, 7.61095, -0.124578, 409164, -330.92, 2707.82, 1.28282e+11,
		1.28282e+11, 258.629, 258.629, 1.17904e-12, -10.4849, 0.253768, 7.86472,
		-0.0991765, 7.16696, 2902.31, 2902.31, 1.33376e+11, 1.33376e+11, 275.96,
		275.96, 8.9308e-05, 8.9308e-05, -7.86477, -7.86472, 2.87394, 2.87394,
		3038.74, 3038.74, 1.37678e+11, 1.37678e+11, 288.117, 288.117,
		9.5024e-05, 9.5024e-05};
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 3);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1500.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 130654281, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 130654281, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.0, 0.01);
	const auto& reg = cluster.getRegion();
	Composition lo = reg.getOrigin();
	Composition hi = reg.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 2);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::Xe] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.2735, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1500.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 0.007707, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.007707, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.0, 0.01);
	const auto& regBis = cluster.getRegion();
	lo = regBis.getOrigin();
	hi = regBis.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::I] = 0;
	comp[Spec::V] = 2;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.386787, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1500.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 2.78553, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 2.78553, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.0, 0.01);
	const auto& regTer = cluster.getRegion();
	lo = regTer.getOrigin();
	hi = regTer.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 3);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::Xe] = 2;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 10);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.44235, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1500.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& regQuad = cluster.getRegion();
	lo = regQuad.getOrigin();
	hi = regQuad.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 3);
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 3);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	// Remove reaction file
	std::remove(reactionFile.c_str());
}

// BOOST_AUTO_TEST_CASE(grouped)
//{
//	// Create the option to create a network
//	xolotl::options::ConfOptions opts;
//	// Create a good parameter file
//	std::string parameterFile = "param.txt";
//	std::ofstream paramFile(parameterFile);
//	paramFile << "netParam=25 0 0 0 0" << std::endl
//			  << "grouping=11 4" << std::endl
//			  << "process=reaction" << std::endl;
//	paramFile.close();
//
//	// Create a fake command line to read the options
//	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
//	opts.readParams(cl.argc, cl.argv);
//
//	std::remove(parameterFile.c_str());
//
//	using NetworkType = NEReactionNetwork;
//	using Spec = NetworkType::Species;
//	using Composition = NetworkType::Composition;
//
//	// Get the boundaries from the options
//	NetworkType::AmountType maxXe = opts.getMaxImpurity();
//	NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
//	// Adapt maxXe
//	int i = 0;
//	while (maxXe + 1 > pow(groupingWidth, i)) {
//		++i;
//	}
//	maxXe = pow(groupingWidth, i) - 1;
//	NetworkType network({maxXe}, {{groupingWidth}}, 1, opts);
//
//	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 16);
//	BOOST_REQUIRE_EQUAL(network.getDOF(), 19);
//
//	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
//	BOOST_REQUIRE_EQUAL(bounds.size(), 16);
//
//	// Get the diagonal fill
//	const auto dof = network.getDOF();
//	NetworkType::SparseFillMap knownDFill;
//	knownDFill[0] = {
//		0, 1, 16, 2, 3, 4, 17, 5, 18, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
//	knownDFill[1] = {1, 0, 16, 4, 17, 5, 18};
//	knownDFill[2] = {2, 0, 8, 9};
//	knownDFill[3] = {3, 0, 11, 12};
//	knownDFill[4] = {4, 0, 17, 14, 1, 16};
//	knownDFill[5] = {5, 0, 1, 16, 18, 6};
//	knownDFill[6] = {6, 0, 5, 18, 15};
//	knownDFill[7] = {7, 0, 8};
//	knownDFill[8] = {8, 0, 7, 2};
//	knownDFill[9] = {9, 0, 2, 10};
//	knownDFill[10] = {10, 0, 9, 11};
//	knownDFill[11] = {11, 0, 10, 3};
//	knownDFill[12] = {12, 0, 3, 13};
//	knownDFill[13] = {13, 0, 12, 14};
//	knownDFill[14] = {14, 0, 13, 4, 17};
//	knownDFill[15] = {15, 0, 6};
//	knownDFill[16] = {16, 0, 1, 4, 17, 5, 18};
//	knownDFill[17] = {17, 0, 4, 14, 1, 16};
//	knownDFill[18] = {18, 0, 1, 16, 5, 6};
//	NetworkType::SparseFillMap dfill;
//	auto nPartials = network.getDiagonalFill(dfill);
//	BOOST_REQUIRE_EQUAL(nPartials, 105);
//	for (NetworkType::IndexType i = 0; i < dof; i++) {
//		auto rowIter = dfill.find(i);
//		if (rowIter != dfill.end()) {
//			const auto& row = rowIter->second;
//			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
//		}
//	}
//
//	// Set temperatures
//	std::vector<double> temperatures = {1000.0};
//	std::vector<double> depths = {1.0};
//	network.setTemperatures(temperatures, depths);
//	NetworkType::IndexType gridId = 0;
//
//	// Check the largest rate
//	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 4.15126e+103, 0.01);
//
//	// Create a concentration vector where every field is at 1.0
//	std::vector<double> concentrations(dof + 1, 1.0);
//	using HostUnmanaged =
//		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
//	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
//	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
//	deep_copy(dConcs, hConcs);
//
//	// Check the total concentrations
//	BOOST_REQUIRE_CLOSE(
//		network.getTotalConcentration(dConcs, Spec::Xe, 1), 25.0, 0.01);
//	BOOST_REQUIRE_CLOSE(
//		network.getTotalAtomConcentration(dConcs, Spec::Xe, 3), 322.0, 0.01);
//	BOOST_REQUIRE_CLOSE(
//		network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1), 16.0555,
//		0.01);
//	BOOST_REQUIRE_CLOSE(
//		network.getTotalTrappedAtomConcentration(dConcs, Spec::Xe, 1), 0.0,
//		0.01);
//
//	// Check the left side rate of 0th cluster
//	BOOST_REQUIRE_CLOSE(
//		network.getLeftSideRate(dConcs, 0, gridId), 0.00010499, 0.01);
//
//	// Create a flux vector where every field is at 0.0
//	std::vector<double> fluxes(dof + 1, 0.0);
//	using HostUnmanaged =
//		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
//	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
//	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
//	deep_copy(dFluxes, hFluxes);
//
//	// Check the fluxes computation
//	std::vector<double> knownFluxes = {-2.07563e+103, -1.76185e-07, -1.8962e-07,
//		-1.13743e-07, 5.18908e+102, -1.5084e-07, 7.36301e-06, 2.45933e-06,
//		-2.38186e-07, -1.60127e-07, -1.39969e-07, -1.25163e-07, -1.04618e-07,
//		-9.71289e-08, -2.07563e+103, 5.1366e-06, -5.49878e-06, -6.22689e+102,
//		-5.81881e-06, 0};
//	network.computeAllFluxes(dConcs, dFluxes, gridId);
//	deep_copy(hFluxes, dFluxes);
//	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
//		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
//	}
//
//	// Check the partials computation
//	std::vector<double> knownPartials = {-0.000121415, -1.9034e-05, 5.44285e-41,
//		-3.44213e-06, -3.98113e-06, 4.15126e+103, -6.22689e+103, -1.99994e-05,
//		4.65987e-41, -5.1366e-06, -3.01432e-06, -3.25251e-06, -3.60226e-06,
//		-3.74222e-06, -3.86739e-06, -4.08575e-06, -4.18288e-06, -4.27373e-06,
//		2.2293e-36, -1.18963e-06, -1.76185e-07, -1.78444e-06, 1.11915e-06,
//		1.67873e-06, 1.53129e-40, -2.29693e-40, -3.44213e-06, -1.8962e-07,
//		3.25251e-06, 6.75237e-17, -3.98113e-06, -1.13743e-07, 3.86739e-06,
//		3.98904e-22, -1.03782e+103, -1.72945e-06, 1.55672e+103, 1.06843e-06,
//		1.44058e-40, -2.16086e-40, -1.24996e-06, -1.5084e-07, 1.18963e-06,
//		1.78444e-06, -1.87494e-06, 2.43742e-49, -5.1366e-06, 7.36301e-06,
//		4.99984e-06, 7.49977e-06, 2.2293e-36, -3.01432e-06, 7.93299e-06,
//		6.13381e-15, -3.25251e-06, -2.38186e-07, 3.01432e-06, 3.6376e-16,
//		-3.60226e-06, -1.60127e-07, 3.44213e-06, 3.95562e-17, -3.74222e-06,
//		-1.39969e-07, 3.60226e-06, 1.26426e-18, -3.86739e-06, -1.25163e-07,
//		3.74222e-06, 4.09399e-19, -4.08575e-06, -1.04618e-07, 3.98113e-06,
//		2.37715e-20, -4.18288e-06, -9.71289e-08, 4.08575e-06, 4.26867e-21,
//		-4.27373e-06, -9.08512e-08, 4.18288e-06, 4.15126e+103, -6.22689e+103,
//		-2.2293e-36, 5.1366e-06, 5.1366e-06, -3.56888e-06, -5.49878e-06,
//		1.42755e-06, -1.34298e-06, -2.01447e-06, 1.83755e-40, -2.75632e-40,
//		-1.86807e+103, -3.29659e-06, 1.24538e+103, -1.28212e-06, 1.72869e-40,
//		-2.59304e-40, -3.74988e-06, -5.81881e-06, -1.42755e-06, -2.14133e-06,
//		1.49995e-06, 2.92491e-49};
//	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
//	network.computeAllPartials(dConcs, vals, gridId);
//	auto hPartials = create_mirror_view(vals);
//	deep_copy(hPartials, vals);
//	int startingIdx = 0;
//	for (NetworkType::IndexType i = 0; i < dof; i++) {
//		auto rowIter = dfill.find(i);
//		if (rowIter != dfill.end()) {
//			const auto& row = rowIter->second;
//			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
//				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
//				auto index = std::distance(row.begin(), iter);
//				XOLOTL_REQUIRE_CLOSE(hPartials[startingIdx + index],
//					knownPartials[startingIdx + j], 0.01);
//			}
//			startingIdx += row.size();
//		}
//	}
//
//	// Check clusters
//	NetworkType::Composition comp = NetworkType::Composition::zero();
//	comp[Spec::Xe] = 1;
//	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
//	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
//	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.0, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
//	BOOST_REQUIRE_CLOSE(
//		cluster.getDiffusionCoefficient(gridId), 3.6298e-7, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1.0, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.0, 0.01);
//	const auto& reg = cluster.getRegion();
//	Composition lo = reg.getOrigin();
//	Composition hi = reg.getUpperLimitPoint();
//	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 1);
//	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 2);
//	auto momId = cluster.getMomentIds();
//	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
//
//	comp[Spec::Xe] = 20;
//	cluster = network.findCluster(comp, plsm::HostMemSpace{});
//	BOOST_REQUIRE_EQUAL(cluster.getId(), 5);
//	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.7961, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 79.0, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
//	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
//	BOOST_REQUIRE_EQUAL(
//		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
//	const auto& regPen = cluster.getRegion();
//	lo = regPen.getOrigin();
//	hi = regPen.getUpperLimitPoint();
//	BOOST_REQUIRE_EQUAL(lo[Spec::Xe], 20);
//	BOOST_REQUIRE_EQUAL(hi[Spec::Xe], 24);
//	momId = cluster.getMomentIds();
//	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
// }

BOOST_AUTO_TEST_CASE(fullyRefined_ReSo)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=5 0 0 5 1" << std::endl
			  << "reactionFilePath=reactionRates.txt" << std::endl
			  << "process=reaction resolution" << std::endl
			  << "tempParam=1500" << std::endl;
	paramFile.close();

	// Create the reaction file
	std::string reactionFile = "reactionRates.txt";
	std::ofstream reacFile(reactionFile);
	reacFile << "0 0 1 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
			 << std::endl
			 << "0 1 0 -0.160789024638881 -12.7387253625273 4.2488E+11 4.23317"
			 << std::endl
			 << "0 2 0 2.87060634573234 -8.37649260139129 9.9337E+11 3.43830674"
			 << std::endl
			 << "1 0 0 9.18245700000039 3.7599999999763 2.9921E+11 1"
			 << std::endl
			 << "1 1 0 3.14296517286615 -1.4324963007117 0 0" << std::endl
			 << "1 2 0 2.26277914822638 -9.53632166324027 6.4829E+11 3.69567683"
			 << std::endl
			 << "1 3 0 4.17605871610343 -5.00715984031393 7.4802E+10 5.34571352"
			 << std::endl
			 << "Reactions" << std::endl
			 << "0 0 1 0.0001" << std::endl
			 << "0 1 0 0.0001" << std::endl
			 << "0 2 0 0.0001" << std::endl;
	reacFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create the network
	using NetworkType = NEReactionNetwork;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxXe = opts.getMaxImpurity();
	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxXe, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{maxXe + 1, maxV + 1, maxI + 1}};
	NetworkType network(maxSpeciesAmounts, subdivRatios, 1, opts);

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 29);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 29);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.547, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0409168, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == true);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 29);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 2);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 20,
		21, 22, 23, 25, 26, 27, 28, 3, 9, 14, 19, 24};
	knownDFill[1] = {1, 0, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 17, 19,
		20, 21, 22, 24, 25, 26, 27, 8, 13, 18, 23, 28};
	knownDFill[2] = {2, 0, 1, 3, 4, 5, 6, 9, 10, 11, 14, 15, 16, 19, 20, 21, 24,
		25, 26, 7, 8, 12, 13, 17, 18, 22, 23, 27, 28};
	knownDFill[3] = {3, 0, 4, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	knownDFill[4] = {4, 0, 5, 1, 3, 2, 6, 9, 11, 12, 10, 13, 15, 16, 17, 18, 20,
		21, 22, 23, 25, 26, 27, 28};
	knownDFill[5] = {5, 0, 6, 1, 4, 2, 3, 9, 10, 11, 14, 15, 16, 19, 20, 21, 7};
	knownDFill[6] = {6, 0, 7, 1, 5, 2, 4, 3, 9, 10, 14, 15, 19, 20, 8, 11, 12};
	knownDFill[7] = {7, 0, 8, 1, 6, 2, 5, 3, 12, 13};
	knownDFill[8] = {8, 0, 1, 7, 2, 6, 3, 13};
	knownDFill[9] = {9, 0, 10, 1, 2, 3, 4, 5, 6, 11, 14, 15};
	knownDFill[10] = {10, 0, 11, 1, 9, 2, 3, 5, 6, 12, 15, 16};
	knownDFill[11] = {11, 0, 12, 1, 10, 2, 9, 3, 6, 4, 5, 13, 16, 17};
	knownDFill[12] = {12, 0, 13, 1, 11, 2, 10, 3, 7, 4, 6, 5, 17, 18};
	knownDFill[13] = {13, 0, 1, 12, 2, 11, 3, 8, 5, 6, 18};
	knownDFill[14] = {14, 0, 15, 1, 2, 3, 9, 5, 6, 16, 19, 20};
	knownDFill[15] = {15, 0, 16, 1, 14, 2, 3, 10, 5, 6, 17, 20, 21};
	knownDFill[16] = {16, 0, 17, 1, 15, 2, 14, 3, 11, 5, 9, 18, 21, 22};
	knownDFill[17] = {17, 0, 18, 1, 16, 2, 15, 3, 12, 5, 10, 6, 9, 22, 23};
	knownDFill[18] = {18, 0, 1, 17, 2, 16, 3, 13, 5, 11, 6, 10, 23};
	knownDFill[19] = {19, 0, 20, 1, 2, 3, 14, 5, 6, 21, 24, 25};
	knownDFill[20] = {20, 0, 21, 1, 19, 2, 3, 15, 5, 6, 22, 25, 26};
	knownDFill[21] = {21, 0, 22, 1, 20, 2, 19, 3, 16, 5, 14, 23, 26, 27};
	knownDFill[22] = {22, 0, 23, 1, 21, 2, 20, 3, 17, 5, 15, 6, 14, 27, 28};
	knownDFill[23] = {23, 0, 1, 22, 2, 21, 3, 18, 5, 16, 6, 15, 28};
	knownDFill[24] = {24, 0, 25, 1, 2, 3, 19, 26};
	knownDFill[25] = {25, 0, 26, 1, 24, 2, 3, 20, 27};
	knownDFill[26] = {26, 0, 27, 1, 25, 2, 24, 3, 21, 5, 19, 28};
	knownDFill[27] = {27, 0, 28, 1, 26, 2, 25, 3, 22, 5, 20, 6, 19};
	knownDFill[28] = {28, 0, 1, 27, 2, 26, 3, 23, 5, 21, 6, 20};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 443);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1500.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 6821697397703317, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::Xe, 1), 26.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::Xe, 3), 60.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Xe, 1), 12.106857,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Xe, 1), 75.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 2704.45, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {2.36776e+06, -1.1276e+11, -1.1276e+11,
		2.57526e+16, 1.44634e+11, -9.12248e+09, -1.28282e+11, -1.33376e+11,
		-1.37678e+11, 6.31097e+15, 385105, -360.232, 3454.38, 3321.87, -390996,
		1.56336e+11, -360.231, 3178.43, 3321.87, -390996, 6.82154e+15, 528872,
		3178.42, 3321.87, -6.31101e+15, -6.82158e+15, 1.28282e+11, 1.33376e+11,
		1.37678e+11, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-2704.45, -8.50534, -2317.66, -6.39485,
		-230.224, -7.31052, -7.61095, -7.86472, -6.93297, -7.31045, -7.61095,
		-7.86472, -6.93296, -7.31045, -7.61095, -7.86472, 409154, -7.31045,
		-7.61095, -7.86472, 409154, -2.91825, -7.61095, -7.86472, 2.48123e-15,
		388034, 388034, 388034, 388034, -1.1276e+11, 2309.16, 2317.67,
		-1.1276e+11, -2.47254, -261.278, -2.70373, -2.87317, -2.47259, -2.65018,
		-2.77463, -2.87338, -2.47259, -2.65018, -2.77476, -2.8739, -2.47259,
		-2.65018, -2.77476, -2.87394, -2.47259, -2.65018, -2.77476, -2.87394,
		0.000781093, 0.000738891, 0.000305525, 1.13101e-05, 7.41541e-09,
		-1.1276e+11, -2317.66, 9.89035, -1.1276e+11, -2707.82, -3160.94,
		-3038.74, -2707.82, -2902.31, -3038.74, -2707.82, -2902.31, -3038.74,
		-2707.82, -2902.31, -3038.74, -2707.82, -2902.31, -3038.74, 4.9183e-09,
		4.61987e-09, 2.36938e-09, 4.13414e-09, 9.24533e-12, 7.06832e-10,
		1.12611e-12, 1.17904e-12, 1.12611e-12, 1.17904e-12, -2.76144e+12,
		6.39485, -1.1276e+11, -1.1276e+11, -1.1276e+11, -1.21883e+11,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 3.18738e+10, -1.21883e+11,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 6.31101e+15, -1.21883e+11,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 6.31101e+15, 3.44525e+10,
		-1.28282e+11, -1.33376e+11, -1.37678e+11, 6.31112e+15, 6.8217e+15,
		120071, 5.17505e-05, 5.34197e-05, -1.1276e+11, 223.83, -28.4037,
		1.1276e+11, 2.19519, -2707.82, -8.93039e-05, 1.44634e+11, 1.00349e-13,
		3.46518e-20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.21883e+11,
		-222.914, -280.736, -258.806, -256.156, 1.1276e+11, -9.12247e+09,
		-258.629, -275.96, -288.117, -258.629, -275.96, -288.117, -258.629,
		-275.96, -288.117, 4.9183e-09, -1.28282e+11, 0.300431, 7.61172, 258.504,
		-26.8386, -330.92, 2707.82, -1.28282e+11, -8.9308e-05, -9.5024e-05,
		-8.9308e-05, -9.5024e-05, -8.9308e-05, -9.5024e-05, 4.61987e-09,
		4.97738e-05, 3.46518e-20, -1.33376e+11, 0.253768, 7.8655, -0.0990872,
		2.77485, 3160.94, 3160.94, -1.33376e+11, 5.17505e-05, 0, -1.37678e+11,
		-7.86472, 2.87394, 2.87394, 3038.74, 3038.74, -1.37678e+11, 5.34197e-05,
		-2.57395e+11, 6.93297, 6.93297, -2.47259, -2707.82, 0, 1.1276e+11,
		-258.629, -8.9308e-05, 1.43843e-10, 6.31112e+15, 0, -1.21883e+11,
		0.377476, 7.31058, -0.177594, 388036, -2902.31, 231.557, 1.21883e+11,
		-9.5024e-05, 2.36938e-09, 17.7205, 0, -1.28282e+11, 0.300508, 7.61151,
		-0.124578, 2.65018, -330.92, 2707.82, 7.62939e-05, 1.28282e+11, 258.629,
		-29.4887, 4.13414e-09, 4.97738e-05, 0, -1.33376e+11, 0.253768, 7.86546,
		-0.0991765, 2.77476, 2902.31, 2902.31, 0, 1.33376e+11, 8.9308e-05,
		8.9308e-05, 1103.84, 5.17505e-05, 0, -1.37678e+11, -7.86472, 2.87394,
		2.87394, 3038.74, 3038.74, 0, 1.37678e+11, 288.117, 288.117,
		5.34197e-05, -6.31123e+15, 6.93297, 6.93297, -2.47259, -2707.82, 0,
		1.1276e+11, -258.629, -8.9308e-05, 1.05064e-12, 6.31112e+15, 0,
		-1.21883e+11, 0.377476, 7.31045, -0.177594, 388036, -2902.31, 0,
		1.21883e+11, -275.96, -9.5024e-05, 9.24533e-12, 1.56335e+11, 0,
		-1.28282e+11, 0.300508, 7.61099, -0.124578, 2.65019, -330.92, 2707.82,
		0, 1.28282e+11, -29.4887, 258.629, 7.06832e-10, 0.00120217, 0,
		-1.33376e+11, 0.253768, 7.86503, -0.0991765, 2.77476, 2902.31, 2902.31,
		0, 1.33376e+11, 275.96, 275.96, 8.9308e-05, 8.9308e-05, 5.17505e-05, 0,
		-1.37678e+11, -7.86472, 2.87394, 2.87394, 3038.74, 3038.74, 0,
		1.37678e+11, 288.117, 288.117, 9.5024e-05, 9.5024e-05, 5.34197e-05,
		-6.31123e+15, 6.93297, 6.93297, -2.47259, -2707.82, 0, 1.1276e+11,
		-258.629, -8.9308e-05, 1.05064e-12, 6.31112e+15, 0, -2.78218e+11,
		0.377476, 7.31045, -0.177594, 388036, -2902.31, 0, 1.21883e+11, -275.96,
		-9.5024e-05, 1.12611e-12, 6.8217e+15, 0, -1.28282e+11, 0.300508,
		7.61095, -0.124578, 409164, -330.92, 2707.82, 0, 1.28282e+11, -29.4887,
		258.629, 1.17904e-12, 120071, 0, -1.33376e+11, 0.253768, 7.86473,
		-0.0991765, 2.77476, 2902.31, 2902.31, 0, 1.33376e+11, 275.96, 275.96,
		8.9308e-05, 8.9308e-05, 5.17505e-05, 0, -1.37678e+11, -7.86472, 2.87394,
		2.87394, 3038.74, 3038.74, 0, 1.37678e+11, 288.117, 288.117, 9.5024e-05,
		9.5024e-05, 5.34197e-05, -6.31112e+15, 6.93297, 6.93297, -2.47259,
		-2707.82, 1.1276e+11, 1.1276e+11, 1.05064e-12, -6.8217e+15, 0.377476,
		7.31045, -0.177594, 388036, -2902.31, 1.21883e+11, 1.21883e+11,
		1.12611e-12, -123124, 0.300508, 7.61095, -0.124578, 409164, -330.92,
		2707.82, 1.28282e+11, 1.28282e+11, 258.629, 258.629, 1.17904e-12,
		-10.4849, 0.253768, 7.86472, -0.0991765, 7.16696, 2902.31, 2902.31,
		1.33376e+11, 1.33376e+11, 275.96, 275.96, 8.9308e-05, 8.9308e-05,
		-7.86477, -7.86472, 2.87394, 2.87394, 3038.74, 3038.74, 1.37678e+11,
		1.37678e+11, 288.117, 288.117, 9.5024e-05, 9.5024e-05};
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

// BOOST_AUTO_TEST_CASE(grouped_ReSo)
//{
//	// Create the option to create a network
//	xolotl::options::ConfOptions opts;
//	// Create a good parameter file
//	std::string parameterFile = "param.txt";
//	std::ofstream paramFile(parameterFile);
//	paramFile << "netParam=25 0 0 0 0" << std::endl
//			  << "grouping=11 4" << std::endl
//			  << "process=reaction resolution" << std::endl;
//	paramFile.close();
//
//	// Create a fake command line to read the options
//	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
//	opts.readParams(cl.argc, cl.argv);
//
//	std::remove(parameterFile.c_str());
//
//	using NetworkType = NEReactionNetwork;
//	using Spec = NetworkType::Species;
//	using Composition = NetworkType::Composition;
//
//	// Get the boundaries from the options
//	NetworkType::AmountType maxXe = opts.getMaxImpurity();
//	NetworkType::AmountType groupingWidth = opts.getGroupingWidthA();
//	// Adapt maxXe
//	int i = 0;
//	while (maxXe + 1 > pow(groupingWidth, i)) {
//		++i;
//	}
//	maxXe = pow(groupingWidth, i) - 1;
//	NetworkType network({maxXe}, {{groupingWidth}}, 1, opts);
//	network.setFissionRate(8.0e-9);
//
//	// Get the diagonal fill
//	const auto dof = network.getDOF();
//	NetworkType::SparseFillMap knownDFill;
//	knownDFill[0] = {
//		0, 1, 16, 2, 3, 4, 17, 5, 18, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
//	knownDFill[1] = {1, 0, 16, 4, 17, 5, 18};
//	knownDFill[2] = {2, 0, 8, 9};
//	knownDFill[3] = {3, 0, 11, 12};
//	knownDFill[4] = {4, 0, 17, 14, 1, 16};
//	knownDFill[5] = {5, 0, 1, 16, 18, 6};
//	knownDFill[6] = {6, 0, 5, 18, 15};
//	knownDFill[7] = {7, 0, 8};
//	knownDFill[8] = {8, 0, 7, 2};
//	knownDFill[9] = {9, 0, 2, 10};
//	knownDFill[10] = {10, 0, 9, 11};
//	knownDFill[11] = {11, 0, 10, 3};
//	knownDFill[12] = {12, 0, 3, 13};
//	knownDFill[13] = {13, 0, 12, 14};
//	knownDFill[14] = {14, 0, 13, 4, 17};
//	knownDFill[15] = {15, 0, 6};
//	knownDFill[16] = {16, 0, 1, 4, 17, 5, 18};
//	knownDFill[17] = {17, 0, 4, 14, 1, 16};
//	knownDFill[18] = {18, 0, 1, 16, 5, 6};
//	NetworkType::SparseFillMap dfill;
//	auto nPartials = network.getDiagonalFill(dfill);
//	BOOST_REQUIRE_EQUAL(nPartials, 105);
//	for (NetworkType::IndexType i = 0; i < dof; i++) {
//		auto rowIter = dfill.find(i);
//		if (rowIter != dfill.end()) {
//			const auto& row = rowIter->second;
//			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
//		}
//	}
//
//	// Set temperatures
//	std::vector<double> temperatures = {1000.0};
//	std::vector<double> depths = {1.0};
//	network.setTemperatures(temperatures, depths);
//	NetworkType::IndexType gridId = 0;
//
//	// Check the largest rate
//	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 8.942685e+107, 0.01);
//
//	// Create a concentration vector where every field is at 1.0
//	std::vector<double> concentrations(dof + 1, 1.0);
//	using HostUnmanaged =
//		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
//	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
//	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
//	deep_copy(dConcs, hConcs);
//
//	// Check the left side rate of 0th cluster
//	BOOST_REQUIRE_CLOSE(
//		network.getLeftSideRate(dConcs, 0, gridId), 2.2618, 0.01);
//
//	// Create a flux vector where every field is at 0.0
//	std::vector<double> fluxes(dof + 1, 0.0);
//	using HostUnmanaged =
//		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
//	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
//	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
//	deep_copy(dFluxes, hFluxes);
//
//	// Check the fluxes computation
//	std::vector<double> knownFluxes = {-4.47134e+107, -0.00379371, -0.00410766,
//		-0.00246059, 1.11784e+107, -0.00318435, 0.158612, 0.0529346,
//		-0.00516169, -0.00346739, -0.00302979, -0.00270841, -0.00226261,
//		-0.00210017, -4.47134e+107, 0.110484, -0.118683, -1.3414e+107,
//		-0.125485, 0};
//	network.computeAllFluxes(dConcs, dFluxes, gridId);
//	deep_copy(hFluxes, dFluxes);
//	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
//		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
//	}
//
//	// Check the partials computation
//	std::vector<double> knownPartials = {-2.61554, -0.409267, 0, -0.0738305,
//		-0.0855092, 8.94269e+107, -1.3414e+108, -0.430117, 0, -0.110482,
//		-0.064144, -0.069715, -0.0773028, -0.0803359, -0.0830468, -0.0877732,
//		-0.0898745, -0.0918394, 0.000168841, -0.0256749, -0.00379539,
//		-0.0383688, 0.0241089, 0.0361633, 4.44896e-05, -6.67345e-05, -0.0744708,
//		-0.00408481, 0.0700659, 0.000297336, -0.0860146, -0.00245026, 0.0833116,
//		0.000242396, -2.23567e+107, -0.0372559, 3.35351e+107, 0.0230163,
//		4.78367e-05, -7.1755e-05, -0.0269713, -0.0032494, 0.025627, 0.0384405,
//		-0.0403234, 4.28086e-05, -0.110824, 0.158615, 0.107707, 0.161561,
//		0.000168841, -0.0653302, 0.170893, 0.000350842, -0.0704167, -0.00513103,
//		0.0649348, 0.000320188, -0.0778975, -0.00344947, 0.0741507, 0.000279413,
//		-0.0808948, -0.00301523, 0.0776001, 0.000264854, -0.0835765,
//		-0.00269627, 0.0806154, 0.000252718, -0.088258, -0.00225369, 0.0857619,
//		0.000233478, -0.0903414, -0.00209236, 0.0880156, 0.000225672,
//		-0.0922907, -0.00195713, 0.0901079, 8.94269e+107, -1.3414e+108,
//		-0.000168841, 0.110653, 0.110653, -0.0770246, -0.118455, 0.030695,
//		-0.0289306, -0.043396, 5.33876e-05, -8.00813e-05, -4.02421e+107,
//		-0.0710155, 2.68281e+107, -0.0276195, 5.7404e-05, -8.6106e-05,
//		-0.0809138, -0.125349, -0.0307524, -0.0461286, 0.0322587, 5.13703e-05};
//	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
//	network.computeAllPartials(dConcs, vals, gridId);
//	auto hPartials = create_mirror_view(vals);
//	deep_copy(hPartials, vals);
//	int startingIdx = 0;
//	for (NetworkType::IndexType i = 0; i < dof; i++) {
//		auto rowIter = dfill.find(i);
//		if (rowIter != dfill.end()) {
//			const auto& row = rowIter->second;
//			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
//				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
//				auto index = std::distance(row.begin(), iter);
//				XOLOTL_REQUIRE_CLOSE(hPartials[startingIdx + index],
//					knownPartials[startingIdx + j], 0.01);
//			}
//			startingIdx += row.size();
//		}
//	}
// }

BOOST_AUTO_TEST_SUITE_END()
