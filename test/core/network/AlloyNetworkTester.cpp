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
	paramFile << "netParam=8 8 0 6 4" << std::endl
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

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 23);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 23);

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
	BOOST_REQUIRE_EQUAL(nPartials, 401);
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
		network.getTotalConcentration(dConcs, Spec::PerfectV, 1), 2.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::PerfectI, 3), 30.0,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::FaultedI, 6),
		1.394446, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::FaultedV, 1),
		0.0, 0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 7043583050, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-1.05369e+09, 2.12557e+09, 3.66059e+09,
		1.06311e+10, -2.18613e+10, -9.32048e+09, -4.80299e+09, 2.45764e+07,
		8.31735e+09, -7.44588e+10, -3.26374e+10, -2.11206e+10, -1.27326e+10,
		-8.14156e+09, -1.24599e+10, -5.60022e+09, -9.33467e+09, -2.83027e+08,
		-1.40243e+08, 1.39461e+09, 2.0646e+09, 1.36196e+09, 3.30628e+09, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-7.04358e+09, -1.71191e+08,
		4.67697e+08, 1.07736e+09, 8712.16, 8360.78, 8062.23, -80436.6, 297817,
		160139, 112625, 3.56182e+09, -5.86224e+09, 3.73301e+09, -1.12494e+08,
		2.29006e+08, 8035.56, 7746.79, -107467, -82617.9, -67456.9, 305852,
		163879, 2.24503e+09, 4.01861e+09, -3.95005e+09, 1.9598e+09, 1.9598e+09,
		3.8455e+09, -1.03758e+08, 7511.29, -160139, -110136, -84619.7, -69058.9,
		-58542, 313364, 1.80502e+09, 1.80502e+09, -806684, 1.35014e+09,
		1.35014e+09, 2.01602e+09, 2.01602e+09, 3.94926e+09, 3.94926e+09,
		-313364, -163879, -112625, -86484, -70549.6, -59784.1, 3.45829e+09,
		3.17469e+09, -3.08055e+10, -4.26622e+09, -2.69838e+09, -2.17416e+09,
		-1.91491e+09, -9.78472e+08, -9.98528e+08, -8.4571e+08, -8.58468e+08,
		-1.25051e+08, -2.22549e+08, -2.99122e+08, 5.86249e+08, 5.27139e+08,
		4.77458e+08, 4.36815e+08, 4.36815e+08, -1.4976e+10, 4.26622e+09,
		4.71828e+07, -2.59113e+09, -2.04353e+09, -8.11585e+08, -8.27656e+08,
		-6.94321e+08, -7.04593e+08, -6.711e+08, -1.55139e+08, -2.1271e+08,
		-6.55896e+08, 5.01398e+08, 4.61942e+08, 4.25993e+08, -1.25727e+10,
		2.69838e+09, 1.68815e+08, 4.21903e+09, 9.63549e+06, -6.98001e+08,
		-7.11424e+08, -5.91742e+08, -6.00343e+08, -6.15285e+07, -6.4187e+08,
		-5.45212e+08, -5.59593e+08, -5.71726e+08, 4.3987e+08, 4.12408e+08,
		-9.34322e+09, 2.17416e+09, 2.17416e+09, 2.59113e+09, 2.59113e+09,
		4.2094e+09, -1.06058e+07, -6.15354e+08, -6.26891e+08, -5.2484e+08,
		-4.67617e+07, -4.61942e+08, -9.54917e+08, -4.89812e+08, -5.0024e+08,
		-5.09349e+08, 3.93108e+08, -4.11979e+09, 5.3732e+09, 1.91491e+09,
		2.04353e+09, 2.04353e+09, 2.52956e+09, 2.52956e+09, 4.22e+09, 4.22e+09,
		-5.52318e+08, -5.62443e+08, -3.93108e+08, -4.12408e+08, -4.25993e+08,
		-8.7363e+08, -4.45954e+08, -4.5394e+08, -9.07009e+10, -3.73293e+09,
		-2.19354e+09, -3.84544e+09, -2.34558e+09, -3.9492e+09, -2.46729e+09,
		-3.48964e+09, -1.21891e+09, -3.56314e+09, -3.63767e+09, -3.71065e+09,
		5.52318e+08, -3.39694e+09, -2.66788e+09, -2.72098e+09, -2.78814e+09,
		-4.04625e+09, -3.0297e+09, -3.07791e+09, -1.19014e+09, -1.69712e+09,
		-4.78842e+10, -1.9597e+09, -8.17784e+08, -2.01594e+09, -9.3708e+08,
		69058.9, -1.03037e+09, 59784.1, -1.10633e+09, -1.94869e+09, 9.20152e+07,
		-1.95243e+09, -1.96997e+09, 5.0024e+08, 4.5394e+08, 6.09466e+09,
		-1.50422e+09, -1.49655e+09, -2.06776e+09, -2.11605e+09, -1.58681e+09,
		-1.61091e+09, -4.11028e+08, -3.26639e+10, -1.34999e+09, -3.53904e+08,
		107467, -4.59425e+08, 84619.7, -5.41617e+08, 70549.6, -6.08099e+08,
		-1.46286e+09, 5.00643e+08, -1.42645e+09, 5.45212e+08, 4.89812e+08,
		4.45954e+08, 3.10083e+08, 1.5747e+09, -1.09239e+09, -1.42212e+09,
		-1.45424e+09, -1.09152e+09, -1.10759e+09, -1.07389e+09, -2.44398e+10,
		289105, -1.2609e+08, 156132, -2.21187e+08, 110136, -2.96109e+08, 86484,
		-3.56846e+08, -1.91491e+09, 5.86249e+08, 5.27139e+08, 4.77458e+08,
		4.36815e+08, 1.30296e+08, 1.88709e+09, 3.37907e+09, -1.09264e+09,
		-1.11671e+09, -8.38734e+08, -8.50787e+08, -8.10701e+08, -8.25523e+08,
		-1.21878e+10, -9.78472e+08, -8.11585e+08, -6.98001e+08, -6.15354e+08,
		-5.52318e+08, 9.69882e+07, 4.04625e+09, -2.06776e+09, -1.42212e+09,
		-1.09264e+09, -305852, -1.24602e+10, -9.98528e+08, -8.27656e+08,
		-7.11424e+08, -6.26891e+08, -5.62443e+08, -4.04625e+09, -2.11605e+09,
		-1.45424e+09, -1.11671e+09, 305852, 305852, -8.67882e+09, -8.4571e+08,
		-6.94321e+08, -5.91742e+08, 4.82093e+07, 3.07791e+09, -1.58681e+09,
		-1.09152e+09, -8.38734e+08, 41175.4, 322739, 199146, 199146, 167749,
		167749, -9.33544e+09, -8.58468e+08, -7.04593e+08, -6.00343e+08,
		-5.2484e+08, -3.07791e+09, -1.61091e+09, -1.10759e+09, -8.50787e+08,
		281564, 281564, 193339, 193339, 156734, 156734, 295039, -8.71822e+09,
		-289105, -297817, -305852, -313364, -7.113e+08, -1.1725e+09,
		7.88345e+08, 6.97695e+08, 5.91463e+08, 7.11862e+08, 6.26893e+08,
		5.24841e+08, -3.93108e+08, 1.43806e+08, 2.69309e+09, 1.22277e+08,
		1.45667e+09, 1.08228e+08, 1.03124e+09, 9.80725e+07, 8.10391e+08,
		-316794, -9.1049e+09, -151778, -156132, -160139, -163879, -7.49688e+08,
		8.88826e+08, 8.11585e+08, 6.94321e+08, 2.57489e+08, 7.11424e+08,
		6.00343e+08, -4.61942e+08, -4.12408e+08, 1.00916e+08, 2.76937e+09,
		9.06005e+07, 1.49655e+09, 8.29874e+07, 1.05392e+09, 7.70413e+07,
		8.2533e+08, 566405, -9.34809e+09, -104562, -107467, -110136, -112625,
		1.0476e+09, 2.07111e+09, 1.68444e+09, 8.93547e+08, 8.27656e+08,
		7.04593e+08, -5.45212e+08, -9.54917e+08, -4.25993e+08, 8.03497e+07,
		2.84941e+09, 7.40857e+07, 1.53052e+09, 6.91443e+07, 1.07389e+09,
		1.15769e+09, 28817.5, 146344, -8.75216e+09, -80436.6, -82617.9,
		-84619.7, -86484, 1.857e+09, 2.11524e+09, 1.70925e+09, -6.55896e+08,
		-5.59593e+08, -4.89812e+08, -8.7363e+08, 6.78568e+07, 2.91739e+09,
		6.35758e+07, 1.56033e+09, 1.48212e+09, 1.42212e+09, 1.09152e+09,
		1.17375e+09, 12109, 175376, 314032, -7.8317e+09, -67456.9, -69058.9,
		-70549.6, -5.71726e+08, -5.0024e+08, -4.45954e+08, 5.9317e+07,
		2.9767e+09, 2.12389e+09, 2.06776e+09, 1.58681e+09, 1.50771e+09,
		1.45424e+09, 1.10759e+09, -8.10701e+08, 5671.81, 331016, 26452.3,
		68864.3, -7.40037e+09, -58542, -59784.1, -5.09349e+08, -4.5394e+08,
		4.00227e+09, 3.94926e+09, 3.0297e+09, 2.16663e+09, 2.11605e+09,
		1.61091e+09, -1.07389e+09, -8.25523e+08, 2605.09, 325344, 15301.8,
		208641, 373467};
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
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::V] = 0;
	comp[Spec::I] = 2;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::I] = 0;
	comp[Spec::PerfectV] = 7;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 7);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 8);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::PerfectV] = 0;
	comp[Spec::FaultedI] = 8;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 8);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 9);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::FaultedI] = 0;
	comp[Spec::PerfectI] = 5;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 5);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.321142, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.2480, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 62637676, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 20736000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.5, 0.01);
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 6);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::PerfectI] = 0;
	comp[Spec::FaultedV] = 7;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 7);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 8);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
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
	paramFile << "netParam=20 20 0 6 4" << std::endl
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

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 39);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 47);

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
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 5048430879.9, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::PerfectV, 1), 14.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::PerfectI, 3), 204.0,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::FaultedI, 6),
		9.377118, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::FaultedV, 1),
		0.0, 0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 90513332027, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-8.30696e+10, -8.23928e+09, -1.21853e+09,
		1.46934e+10, 3.27115e+10, -6.71954e+10, -2.67127e+10, -8.78605e+09,
		1.18829e+10, -2.71057e+11, -1.28822e+11, -8.4916e+10, -6.18613e+10,
		-5.13937e+09, -5.24114e+09, -6.64034e+09, -1.99979e+10, 6.17896e+08,
		-7.61419e+08, -3.77277e+09, -1.469e+10, 4.76942e+09, 5.29665e+09,
		7.34094e+09, 8.66039e+09, 1.03921e+10, 9.57801e+09, -6.34134e+09,
		-6.73424e+09, -3.65352e+09, -3.08849e+09, -5.73864e+10, -4.58578e+10,
		-3.53026e+10, -3.06594e+10, -2.51269e+10, -2.12549e+10, -1.77509e+10,
		-1.29565e+10, -6.11325e+09, -2.01642e+09, 1.7262e+09, 4.43767e+09,
		-6.19554e+09, -1.01101e+10, -4.96386e+08, -3.78849e+09, -1.14564e+10,
		-6.9672e+09, -2.10452e+10, -2.21237e+10, -1.31182e+10, -1.53363e+10, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-9.74642e+10, -9.33235e+08,
		-5.33884e+09, 0, -2.28264e+09, 3.42396e+09, -2.60654e+09, -2.34454e+09,
		-1.90016e+09, -4.26622e+09, -2.69838e+09, -2.17416e+09, -1.91491e+09,
		-5.51498e+08, -3.61306e+09, 1.69076e+08, -5.08719e+09, 9.35284e+08,
		-1.18095e+09, -8.4571e+08, -4.50128e+09, -2.38419e-07, -4.70771e+09, 0,
		-9.6263e+08, -1.25051e+08, -2.22549e+08, -2.99122e+08, -9.56868e+08,
		-9.78472e+08, -9.98528e+08, -1.01736e+09, -2.49515e+09, -2.42863e+09,
		-2.38693e+09, -2.36064e+09, -2.25983e+09, -2.25517e+09, -2.25593e+09,
		-1.48818e+09, -2.27282e+09, -5.91372e+08, -6.29027e+08, -8.58468e+08,
		-8.7019e+08, 3.75158e+08, 6.07073e+08, 0, -3.1321e+08, 9.4762e+08,
		4.36815e+08, 4.0326e+08, -1.42292e+10, -9.33235e+08, -7.75244e+08,
		-4.33586e+08, -3.06128e+08, -1.71191e+08, 4.67697e+08, 1.07736e+09,
		2.5232e+09, 8712.16, 8360.78, 8062.23, 7802.01, -6.67579e+08,
		-5.89157e+08, -5.2671e+08, -4.75189e+08, -3.99255e+08, -3.7042e+08,
		-3.45842e+08, -3.2463e+08, 139797, -215512, -75714.7, 0, 297817, 160139,
		112625, 88238.6, -7.39323e+09, 9.32495e+08, -3.82107e+09, 3.71577e+08,
		9.29126e+08, 1.12363e+09, 1.44492e+09, 2.15666e+09, 1.37428e+09,
		1.36209e+09, 9.85582e+08, 7.26965e+08, 6.06776e+08, 1.55856e+09,
		-1.803e+09, -1.42854e+09, -1.00686e+09, -5.57166e+08, -5590.15,
		-2778.61, 0, 2758.5, 5502.88, 32532.7, 444131, -352238, 14225.1,
		-2.45361e+09, 2.60598e+09, 6.78053e+09, 3.69282e+09, -1.52763e+09,
		2.62203e+08, 1.45328e+09, 1.14258e+09, 8.47582e+08, 5.58058e+08,
		3.74747e+08, 1.51732e+09, 2.24716e+09, -2.16545e+08, 1.29868e+09,
		1.54204e+09, 1.44027e+09, 1.53836e+09, 1.29888e+09, 1.4613e+09,
		1.52737e+09, 150004, 77810, 27092.8, -19226.1, -64380.8, -68658.7,
		180801, 4.19538e+08, 4.6796e+08, 5.21438e+08, 5.86287e+08, -975132,
		1.30802e+09, 1.22658e+10, 7.64003e+09, 1.30802e+09, 3.06128e+08,
		1.44857e+09, 3.06139e+09, 6.43009e+08, 1.33873e+09, 3.31242e+09,
		2.30832e+09, 1.53433e+09, 9.19167e+08, 2.00847e+08, 1.42355e+10,
		4.92692e+09, 1.68497e+09, 0, -382527, -198266, -135464, -103565, -84184,
		-71125.3, 3.36026e+08, 1.3759e+09, 3.69365e+08, 1.34485e+09,
		4.07156e+08, 1.34443e+09, 4.4297e+08, 1.36214e+09, -8.37132e+10,
		1.07343e+09, 4.71828e+07, -7.75244e+08, -4.4149e+09, 0, -9.41595e+08,
		1.88319e+09, -2.04879e+09, -1.20208e+09, -2.59113e+09, -2.04353e+09,
		-1.77071e+09, -4.15714e+08, -2.77118e+09, 1.53998e+08, -4.40444e+09,
		6.07073e+08, -9.73158e+08, -6.94321e+08, -3.69121e+09, 2.38419e-07,
		-3.85745e+09, 0, -7.88473e+08, -6.711e+08, -1.55139e+08, -2.1271e+08,
		-2.60766e+08, -7.94246e+08, -8.11585e+08, -8.27656e+08, -8.42722e+08,
		-2.26273e+09, -2.17238e+09, -2.11451e+09, -2.07552e+09, -2.03049e+09,
		-1.9488e+09, -1.93987e+09, -1.93479e+09, -4.47824e+08, -4.77368e+08,
		-7.04593e+08, -7.14032e+08, 3.68759e+08, 9.1061e+08, 6.86424e+08,
		3.27432e+08, 0, -2.98936e+08, 3.46086e+08, -4.7095e+10, -1.20855e+09,
		9.71434e+08, -4.33586e+08, -4.98038e+08, 9.96075e+08, 5.33924e+08,
		-3.42433e+07, 1.28083e+07, 7.55244e+07, 1.60846e+08, -1.68695e+08,
		-2.20416e+09, 5.72063e+08, -2.6355e+09, 0, -5.42666e+08, -1.61246e+09,
		-8.06229e+08, -2.10011e+09, 1.19209e-07, -4.28665e+08, -2.45976e+07,
		-4.79176e+07, -6.91395e+07, -8.83885e+07, -1.05891e+08, -4.97029e+08,
		-4.49453e+08, -4.62112e+08, -4.70265e+08, -4.77871e+08, -1.68759e+08,
		4.98227e+08, 2.37476e+09, -1.83176e+08, -1.96943e+08, 2.62868e+08,
		8.58094e+08, 2.93417e+08, 3.01941e+08, 2.94008e+08, 2.83642e+08,
		2.73052e+08, 2.53287e+08, 2.87089e+08, 2.80928e+08, -2.49331e+10,
		-1.27862e+09, 2.8443e+09, -3.06128e+08, -1.07611e+09, 1.20877e+09,
		-4.14717e+07, -2.1819e+07, 1.95727e+06, 3.0746e+07, -3.26549e+08,
		-2.16561e+09, -3.60934e+08, -1.90161e+09, -1.19209e-07, -3.90648e+08,
		-1.18817e+09, -5.94085e+08, -3.02831e+08, -1.17579e+07, -2.34548e+07,
		-3.45656e+07, -4.5025e+07, -5.48483e+07, -2.6712e+08, -3.35642e+08,
		-3.44172e+08, 2.0802e+08, 4.46401e+09, 2.47648e+09, 1.81749e+09,
		2.2323e+08, 2.21647e+08, 2.1771e+08, 2.12986e+08, 1.13901e+09,
		1.13901e+09, 1.11341e+09, 1.11341e+09, -5.8592e+09, 1.30802e+09,
		2.77603e+09, 1.14244e+09, 1.14244e+09, 1.83676e+09, -2.55659e+08,
		-1.43016e+09, -1.19209e-07, -1.82122e+09, -3.03537e+08, -3.11278e+08,
		-1.87359e+08, -1.94826e+08, -2.00145e+08, -2.04407e+08, -2.0802e+08,
		-2.11186e+08, -2.62495e+08, -2.68913e+08, 4.67163e+09, 4.67163e+09,
		2.53915e+09, 2.53915e+09, 1.82262e+09, 1.82262e+09, 1.46801e+09,
		1.03987e+09, 1.03987e+09, 9.75486e+08, 9.75486e+08, 9.37275e+08,
		9.37275e+08, 9.19167e+08, 9.19167e+08, -2.92376e+11, -3.48964e+09,
		-1.21891e+09, -3.73293e+09, -2.19354e+09, -2.19213e+10, 0, -2.37259e+10,
		0, -3.56314e+09, -3.82844e+09, -3.6599e+09, -4.86642e+09, -1.98837e+10,
		2.5018e+08, -3.03537e+08, -2.2071e+10, 1.44419e+08, -2.66788e+09,
		-2.72098e+09, -2.78814e+09, -3.39694e+09, -4.93931e+09, -3.0297e+09,
		-1.61809e+10, 0, -1.69656e+10, 0, -3.47348e+09, -1.19014e+09,
		-1.69712e+09, -2.34558e+09, -2.46729e+09, -3.84544e+09, -3.9492e+09,
		-4.04625e+09, -4.13777e+09, -3.63767e+09, -3.71065e+09, -3.69009e+09,
		-3.76051e+09, -2.95952e+09, -3.52185e+09, -4.01903e+09, -4.46401e+09,
		-3.46293e+09, -3.89965e+09, -4.30083e+09, -4.67163e+09, -3.5347e+09,
		-3.07791e+09, -3.12225e+09, -1.4963e+11, -1.94869e+09, 9.20152e+07,
		-1.9597e+09, -8.17784e+08, -1.14203e+10, 0, -9.85384e+09, 4.92692e+09,
		-1.95243e+09, -1.99045e+09, -1.64579e+09, -2.4983e+09, -9.50089e+09,
		2.17335e+08, 0, -1.10238e+10, 7.15071e+08, 6.09466e+09, -1.50422e+09,
		-1.49655e+09, -1.56524e+09, -2.56008e+09, -1.58681e+09, -8.45026e+09, 0,
		-8.8426e+09, 0, -1.8087e+09, -4.11028e+08, -9.3708e+08, -1.03037e+09,
		-1.10633e+09, -2.01594e+09, -2.06769e+09, -2.11599e+09, -2.16157e+09,
		-1.96997e+09, -1.99343e+09, -2.01959e+09, -1.96111e+09, -2.01976e+09,
		-1.17713e+09, -1.66864e+09, -2.10565e+09, -2.85431e+09, -1.49109e+09,
		-1.88019e+09, -1.7186e+09, -1.61091e+09, -1.63308e+09, 3.15873e+08,
		-1.02112e+11, -1.46286e+09, 5.00643e+08, -1.34999e+09, -3.53904e+08,
		-7.82931e+09, 0, -5.0549e+09, 5.0549e+09, -1.42645e+09, -1.37178e+09,
		-9.60004e+08, -1.4546e+09, -6.01282e+09, 1.90759e+08, 3.03537e+08,
		-7.47595e+09, 9.74829e+08, 3.10083e+08, 1.5747e+09, -1.09239e+09,
		-1.74915e+09, -1.09152e+09, -5.80176e+09, 0, -6.06332e+09, 0,
		-1.23945e+09, -4.59425e+08, -5.41617e+08, -6.08099e+08, -6.63462e+08,
		-1.38754e+09, -1.42203e+09, -1.45417e+09, -1.48445e+09, -1.41719e+09,
		-1.41961e+09, -1.42795e+09, -1.43955e+09, -1.38755e+09, -1.40396e+09,
		-6.02819e+08, -1.05211e+09, -1.81848e+09, -2.15006e+09, -1.01821e+09,
		-1.07098e+09, -1.10759e+09, -1.12237e+09, 9.88434e+08, 6.31747e+08,
		-7.90189e+10, -1.91491e+09, -1.03833e+09, -1.2609e+08, -6.00104e+09,
		-128165, -4.50429e+08, -2.57854e+09, 3.8678e+09, -1.18446e+09,
		-1.14254e+09, -7.48284e+08, -4.27532e+09, 1.69076e+08, 6.07073e+08,
		-5.82693e+09, 9.35284e+08, 1.30296e+08, 1.88709e+09, 3.37907e+09,
		-6.65662e+08, -1.33727e+09, -8.38734e+08, -4.45178e+09, 2.38419e-07,
		-4.64795e+09, 0, -9.49679e+08, -2.21187e+08, -2.96109e+08, -3.56846e+08,
		-4.07369e+08, -1.06663e+09, -1.09253e+09, -1.11662e+09, -1.13929e+09,
		-1.15141e+09, -1.13866e+09, -1.13543e+09, -1.13742e+09, -1.07392e+09,
		-1.08154e+09, -1.09166e+09, -3.3101e+08, -1.1248e+09, -7.09552e+08,
		-7.51034e+08, -8.50787e+08, -8.61872e+08, 0, -3.1321e+08, 9.4762e+08,
		-1.68497e+10, -9.78472e+08, 6.14076e+09, -3.30937e+09, 7.13953e+07,
		6.49881e+07, 2.67756e+09, -4.44646e+08, -3.26549e+08, -2.55659e+08,
		9.69882e+07, 4.04625e+09, 9.38122e+07, 2.16157e+09, -2.98798e+09,
		-1.09264e+09, -305852, -160139, -110136, -84619.7, -69058.9, -58542,
		8.17031e+08, 1.38389e+09, -1.18577e+09, -5.02524e+08, 1.21532e+08,
		9.67916e+08, 6.88339e+07, -3.46495e+08, -3.08994e+08, -2.93417e+08,
		-2.79498e+08, -2.6698e+08, 4.38199e+08, -1.24407e+10, 2.94205e+07,
		3.26166e+09, 7.2696e+09, -4.00209e+09, 5.86142e+07, -3.89505e+08,
		6.70107e+08, -5.05308e+08, -2.86031e+08, -7.21818e+07, -3.58783e+07, 0,
		3.56187e+07, 64471.4, -77976.8, -59673.5, -33402.7, -40417.5, -139634,
		-110634, 3.76432e+08, 3.29361e+08, 1.66989e+08, -7.91566e+07,
		-8.30393e+08, -9.46852e+08, -8.88623e+08, -9.14394e+08, -4.11855e+08,
		-3.27432e+08, -2.50568e+08, -2.98936e+08, 152603, 98664.2, -1.30106e+10,
		-2.20337e+08, 3.8215e+09, 2.54799e+09, -7.46963e+08, -5.271e+08,
		-3.80322e+08, -4.24951e+08, 1.9369e+09, 1.00471e+09, 3.4983e+08,
		-2.48254e+08, 412768, 326260, -16770.4, 100577, 119424, 111542, 91958.2,
		37823.4, -8.06219e+08, -7.08071e+08, -6.33451e+08, -5.74676e+08,
		-4.8774e+08, -4.54596e+08, -4.2627e+08, -4.01759e+08, -6.50521e+08,
		-7.58042e+08, -7.24851e+08, -7.58096e+08, 12190.3, -1.99997e+10,
		-1.18095e+09, -9.73158e+08, -5.42666e+08, -3.90648e+08, -3.11278e+08,
		-4.93931e+09, -2.56008e+09, -1.74915e+09, -1.33727e+09, 64082.7, 128165,
		192248, 869834, 663484, 1.10248e+06, 381567, 130493, 0, -81217.6,
		-8.32422e+08, -7.30459e+08, -6.52987e+08, -5.92003e+08, -5.0187e+08,
		-4.67532e+08, -4.38199e+08, -4.12828e+08, -3.71077e+08, -3.53671e+08,
		-3.38077e+08, -3.24021e+08, -8.67926e+09, -8.4571e+08, 4.86682e+09,
		-2.38462e+09, 4.39208e+07, 4.20021e+08, 2.13599e+09, -3.74753e+08,
		4.82093e+07, 3.07791e+09, 4.62756e+07, 1.63308e+09, -2.25188e+09,
		-8.38734e+08, 41175.4, 51677.2, 66308.8, 89801.6, 135527, 268804,
		6.65306e+08, 1.64726e+09, -5.09837e+08, 0, 7.73807e+08, 1.07725e+09,
		3.42042e+08, 3.42042e+08, -7.04624e+09, 1.6515e+07, 3.00059e+09,
		5.54693e+09, -2.79806e+09, 3.32491e+07, -3.98004e+08, 5.3498e+08,
		-3.13868e+07, -1.56934e+07, 0, 1.56934e+07, 59388.4, -67668, -26299.8,
		43427.6, 55617.8, -122419, -96768.4, 2.93249e+08, 2.55645e+08,
		1.25161e+08, -7.13624e+07, -4.45816e+08, -3.44964e+08, -1.93396e+08,
		138336, 90181.5, -8.96223e+09, -1.84091e+08, 3.34431e+09, 1.84448e+09,
		-6.13796e+08, -4.20021e+08, -3.56451e+08, 1.37332e+09, 7.15443e+08,
		2.4789e+08, -1.819e+08, 367039, 289487, -8750.75, 91321.6, 106865,
		99293.9, 81532.2, 33927.1, -6.56377e+08, -5.73059e+08, -5.09837e+08,
		-4.60149e+08, -3.86903e+08, -3.59084e+08, -3.35368e+08, -3.14899e+08,
		-3.37585e+08, -2.67367e+08, -1.52937e+08, 11073.7, -1.46915e+10,
		-9.6263e+08, -7.88473e+08, -4.28665e+08, -3.02831e+08, -3.47348e+09,
		-1.8087e+09, -1.23945e+09, -9.49679e+08, 57133.4, 114267, 171400,
		749158, 567893, 946015, 328713, 112698, 0, -70374.3, -6.70604e+08,
		-5.85306e+08, -5.2059e+08, -4.69734e+08, -3.94774e+08, -3.66307e+08,
		-3.42042e+08, -3.21099e+08, -2.86748e+08, -2.72474e+08, -2.59716e+08,
		-2.48241e+08, -1.1765e+10, -7.113e+08, -289105, -1.69771e+06, 0,
		-1.83746e+06, 0, -382527, -1.1725e+09, -3.26539e+08, 2.33434e+09,
		-3.44537e+08, 1.88182e+09, -2.96294e+08, -2.34988e+08, 1.81569e+09,
		-1.80318e+08, 1.41395e+09, -1.49684e+08, -1.87359e+08, 1.43806e+08,
		2.69309e+09, 1.22277e+08, 1.45667e+09, 1.08228e+08, 1.03124e+09,
		9.80725e+07, 8.10391e+08, 6.97695e+08, 7.88345e+08, -1.34006e+09,
		5.30902e+08, 1.26945e+09, 1.91515e+09, -9.64981e+08, 4.19492e+08,
		1.01109e+09, 1.54969e+09, 3.24021e+08, 3.77267e+08, 5.91463e+08,
		2.48241e+08, -316794, -297817, -305852, -313364, -320452, 7.11862e+08,
		6.52281e+08, 6.2658e+08, 5.71586e+08, 5.24555e+08, 4.73195e+08,
		-1.22791e+10, -7.49688e+08, -151778, -884453, 0, -953918, 0, -198266,
		8.88826e+08, 8.11425e+08, 6.94174e+08, 5.59227e+08, 2.52829e+09,
		-4.14377e+08, 2.049e+09, -3.57547e+08, 4.32262e+08, 1.91355e+09,
		-2.03386e+08, 1.49534e+09, -1.69251e+08, -1.94826e+08, 1.00916e+08,
		2.76937e+09, 9.06005e+07, 1.49655e+09, 8.29874e+07, 1.05392e+09,
		7.70413e+07, 8.2533e+08, -1.49826e+09, -3.73434e+08, 1.34979e+09,
		2.03204e+09, -1.03697e+09, -2.59194e+08, 1.05276e+09, 1.61415e+09,
		3.38077e+08, 3.86021e+08, 2.59716e+08, 566405, -156132, -160139,
		-163879, -167404, 2.57489e+08, 7.0738e+08, -2.0288e+08, 7.1126e+08,
		6.37515e+08, 6.00193e+08, 5.31488e+08, -1.26151e+10, 1.0476e+09,
		2.071e+09, 1.68434e+09, -104562, -606344, 0, -652465, 0, -135464,
		8.93547e+08, 8.27543e+08, 7.0449e+08, 1.44952e+09, 2.7646e+09,
		-5.09194e+08, 2.25337e+09, -4.40826e+08, 1.10245e+09, 2.02395e+09,
		-2.31946e+08, 1.58738e+09, -1.93513e+08, -2.00145e+08, 8.03497e+07,
		2.84941e+09, 7.40857e+07, 1.53052e+09, 6.91443e+07, 1.07374e+09,
		1.15769e+09, -2.17042e+09, -4.25993e+08, 5.98255e+08, 2.17278e+09,
		-1.11671e+09, -2.83642e+08, 4.49762e+08, 1.68877e+09, 3.53671e+08,
		3.98197e+08, 2.72474e+08, 28817.5, 146344, -107467, -110136, -112625,
		-114968, 7.87027e+08, -2.1771e+08, -2.08492e+08, 7.23878e+08,
		6.08143e+08, -1.20941e+10, 1.857e+09, 2.11515e+09, 1.70918e+09,
		-80436.6, -464758, -2.91038e-11, -499240, 0, -103565, 9.00858e+08,
		8.42634e+08, 7.13951e+08, 2.34202e+09, 3.06093e+09, -6.43009e+08,
		2.5105e+09, -5.58507e+08, 1.77408e+09, 2.15049e+09, -2.67375e+08,
		1.69312e+09, -2.23647e+08, -2.04407e+08, 6.78568e+07, 2.91726e+09,
		6.35758e+07, 1.5602e+09, 1.48212e+09, 1.42203e+09, 1.09144e+09,
		1.17375e+09, -1.94563e+09, -4.89812e+08, 2.12223e+08, 1.57239e+09,
		-1.20826e+09, -3.10231e+08, 4.71356e+08, 1.16026e+09, 3.71077e+08,
		4.13065e+08, 2.86748e+08, 12109, 175376, 314032, -82617.9, -84619.7,
		-86484, -88238.6, -2.32991e+08, -2.22463e+08, -2.12986e+08,
		-1.16411e+10, 1.88755e+09, 2.15665e+09, 1.732e+09, -378574, 0, -406088,
		0, -84184, -1.62122e+09, 3.44596e+09, -8.41692e+08, 2.84571e+09,
		-7.33416e+08, -1.31588e+09, 2.29738e+09, -3.11991e+08, 1.8162e+09,
		-2.61638e+08, 4.30611e+08, 3.90648e+08, 3.02831e+08, -2.0802e+08,
		5.9317e+07, 2.97659e+09, 2.12389e+09, 2.06769e+09, 1.58674e+09,
		1.50771e+09, 1.45417e+09, 1.10753e+09, 1.19054e+09, -5.71726e+08,
		7.15262e+08, 1.72575e+09, 2.14477e+09, -3.40269e+08, 4.97355e+08,
		1.22738e+09, 1.87171e+09, 5671.81, 195410, 26452.3, 68864.3, -67456.9,
		-69058.9, -70549.6, -71951.7, -2.49306e+08, -2.37263e+08, -2.26491e+08,
		-2.16796e+08, -1.12377e+10, -1.96802e+09, 5.17164e+09, -3.55843e+09,
		4.19497e+09, -2.79279e+09, -320414, -1.45519e-11, -343293, 0, -71125.3,
		0, -7.50316e+08, 2.47038e+09, -3.69317e+08, 1.96153e+09, -3.10499e+08,
		-2.6712e+08, -2.11186e+08, 4.00227e+09, 3.9492e+09, 3.02965e+09,
		2.16663e+09, 2.11599e+09, 1.61086e+09, 1.53299e+09, 1.48445e+09,
		1.12232e+09, -2.9161e+09, 1.38604e+09, 1.92165e+09, 2.80362e+09,
		-1.44519e+09, 5.28557e+08, 1.3057e+09, 1.9841e+09, 4.12828e+08,
		4.51103e+08, 3.21099e+08, 2605.09, 210516, 15301.8, 80782.5, 223772,
		-58542, -59784.1, -60951.6, -2.53287e+08, -2.41005e+08, -2.3002e+08,
		-2.20134e+08, -1.43292e+10, -9.56868e+08, 3.73301e+09, -1.12494e+08,
		232749, -261618, -175570, -67456.9, 64082.7, -7.94246e+08, -4.49453e+08,
		2.29006e+08, 1.17824e+09, 2.63096e+09, 8035.56, 305852, 7746.79, 163879,
		7501.03, 114968, -6.83498e+08, -6.02874e+08, -5.41353e+08, -4.92737e+08,
		-4.13708e+08, -3.83689e+08, -3.58104e+08, -3.36026e+08, -1.43451e+10,
		-9.78472e+08, 1.9598e+09, -1.07957e+08, 354018, -271407, -231405,
		-84619.7, 6655.83, 69623.5, -8.11585e+08, -4.62112e+08, -1.03758e+08,
		3.8455e+09, 3.82907e+08, 7.12386e+08, 7511.29, 313364, 7265.33, 167404,
		-6.98001e+08, -6.15354e+08, -5.52318e+08, -5.02524e+08, -4.26974e+08,
		-3.95871e+08, -3.69365e+08, -1.4286e+10, 2.45976e+09, 1.35014e+09,
		-1.04102e+08, 530909, -271170, -340769, -112625, 6467.66, 80879.9,
		132464, -8.27656e+08, -4.70265e+08, -9.69882e+07, 3.94926e+09,
		-1.00029e+08, 2.01602e+09, 2.05799e+09, 7087.95, 320452, -7.11424e+08,
		-6.26891e+08, -5.62443e+08, -5.11553e+08, -4.36036e+08, -4.07156e+08,
		-1.41585e+10, 1.51227e+09, 1.97186e+09, 7.5367e+08, -1.00742e+08,
		806367, -611592, -659993, -167404, 6300.64, 97664.7, 155192, 68658.7,
		-137317, -129610, -4.77871e+08, -9.15219e+07, 4.04625e+09, -9.38122e+07,
		2.06776e+09, -9.68557e+07, 1.38765e+09, -7.23993e+08, -6.37682e+08,
		-5.71905e+08, -5.19983e+08, -4.4297e+08, -7.45434e+10, -3.23917e+08,
		1.68815e+08, -6.67579e+08, -3.78758e+09, 0, 1.45491e+09, 9.63549e+06,
		-1.85682e+09, -1.98152e+08, 3.44172e+08, -3.79821e+08, -2.86031e+08,
		-2.21876e+09, 1.42858e+08, -1.9624e+09, -1.67854e+09, -3.28636e+08,
		-1.73461e+09, 1.02964e+09, 6.54865e+08, 3.1321e+08, 0, -4.0311e+09, 0,
		-8.32422e+08, -5.91742e+08, -3.14262e+09, -2.38419e-07, -3.28188e+09, 0,
		-6.70604e+08, -6.15285e+07, -2.01822e+09, -6.4187e+08, -1.94304e+09,
		-1.59706e+08, -1.89252e+09, -2.31457e+08, -1.83104e+09, 3.21654e+08,
		-1.81231e+09, -6.83498e+08, -6.98001e+08, -7.11424e+08, -7.23993e+08,
		-1.72352e+09, -3.55198e+08, -6.00343e+08, -6.08247e+08, -6.44912e+10,
		-7.31935e+08, 2.64736e+08, -5.89157e+08, -2.66566e+09, 1.33283e+09,
		-4.31851e+07, 9.7461e+07, -1.7254e+09, -6.02135e+08, 3.35642e+08,
		-2.91249e+08, 0, -2.10647e+09, 6.82044e+08, -1.06058e+07, 1.75131e+09,
		-1.61612e+09, -2.6886e+08, 3.46495e+08, 9.82297e+08, 6.2642e+08,
		2.98936e+08, -3.54036e+09, 0, -7.30459e+08, -2.74543e+09, 0,
		-2.86529e+09, 0, -5.85306e+08, -4.67617e+07, -1.8337e+09, -8.85077e+07,
		-1.76984e+09, -1.56194e+08, -1.69288e+09, -1.83916e+08, -1.66855e+09,
		-2.08532e+08, -1.65016e+09, -6.02874e+08, -6.15354e+08, -6.26891e+08,
		-6.37682e+08, -3.12152e+08, -5.2484e+08, -5.31639e+08, 3.25531e+08,
		-5.66964e+10, 2.54945e+09, 3.42668e+08, -5.2671e+08, -1.79118e+09,
		1.79118e+09, -4.96969e+08, 1.69627e+08, -1.63338e+09, -7.04825e+07,
		3.26549e+08, -2.25768e+08, 2.86031e+08, -2.08712e+09, 9.25357e+08,
		-2.24086e+07, 1.99319e+09, 5.60325e+07, 1.74115e+08, 9.3963e+08,
		5.97873e+08, -3.16725e+09, 0, -6.52987e+08, -2.44397e+09, -1.19209e-07,
		-2.54919e+09, 0, -5.2059e+08, -3.68048e+07, -1.68855e+09, -1.00241e+08,
		-1.59394e+09, -5.63399e+08, -1.56432e+09, -1.50054e+08, 2.959e+08,
		-1.71117e+08, 2.82824e+08, -5.41353e+08, -5.52318e+08, -5.62443e+08,
		-5.71905e+08, -2.44918e+08, -2.62913e+08, -4.73484e+08, 3.17525e+08,
		3.08992e+08, -5.08211e+10, -1.1675e+09, 8.96454e+08, -4.75189e+08,
		-1.08473e+09, 1.62709e+09, 2.30028e+08, -1.57077e+09, -2.97645e+07,
		2.6712e+08, -1.43247e+08, 5.72063e+08, -2.12685e+09, 8.83904e+08,
		-2.96523e+07, 2.19756e+09, 3.00402e+07, 3.50399e+08, 1.12486e+08,
		-3.15625e+08, -1.9353e+08, 3.08994e+08, 8.96809e+08, -2.87338e+09, 0,
		-5.92003e+08, -2.20699e+09, 1.19209e-07, -2.30075e+09, -1.19209e-07,
		-4.69734e+08, -5.75384e+07, -1.52244e+09, -8.24779e+07, 3.03028e+08,
		-1.04848e+08, 2.90282e+08, -5.28257e+08, 2.78263e+08, -4.92737e+08,
		-5.02524e+08, -5.11553e+08, -5.19983e+08, -2.10101e+08, -2.25768e+08,
		3.01452e+08, 2.94215e+08, -4.06163e+10, -1.23376e+09, 1.04266e+09,
		-3.99255e+08, -9.94515e+08, -7.89247e+08, 4.30631e+09, -3.72021e+07,
		2.58011e+08, -7.56069e+07, 2.68913e+08, -1.74056e+08, 907726, 6.253e+08,
		5.03145e+07, -4.57905e+07, 1.15219e+08, -5.74451e+08, 2.79498e+08,
		-2.30665e+09, -1.84532e+09, -2.4387e+09, 1.19209e-07, -5.0187e+08,
		-1.11454e+09, -1.11454e+09, -1.93452e+09, 0, -3.94774e+08, -2.06903e+07,
		2.81251e+08, -4.05646e+07, 2.75331e+08, -5.88561e+07, 2.66896e+08,
		-9.09629e+07, 2.49306e+08, -1.05082e+08, 2.41005e+08, -4.13708e+08,
		-4.26974e+08, -4.36036e+08, -4.4297e+08, -1.61862e+08, 2.74174e+08,
		-3.75261e+10, -1.25045e+09, 1.10957e+09, -3.7042e+08, -1.02751e+09,
		-6.67713e+08, 2.46072e+09, -7.57007e+06, 2.52275e+08, -5.07532e+07,
		2.62495e+08, -5.51818e+08, -3.91126e+07, 4.34351e+09, 3.24198e+07,
		1.84991e+09, 8.35461e+07, -1.28891e+08, 2.6698e+08, -2.15194e+09,
		-2.58232e+09, -2.27298e+09, 1.19209e-07, -4.67532e+08, -6.89927e+08,
		-1.03489e+09, -1.79542e+09, 1.19209e-07, -3.66307e+08, -1.76618e+07,
		2.63589e+08, -3.4814e+07, 2.59194e+08, -6.54714e+07, 2.4476e+08,
		-7.90612e+07, 2.37263e+08, -9.16344e+07, 2.3002e+08, -3.83689e+08,
		-3.95871e+08, -4.07156e+08, 4.64556e+08, -4.1353e+08, -3.35202e+10,
		-1.26231e+09, 2.92263e+09, -3.45842e+08, -1.04945e+09, 1.27025e+09,
		1.8339e+09, 1.93176e+07, 2.45102e+08, -3.02287e+07, 2.55659e+08,
		-1.13706e+08, -4.03264e+07, 4.38262e+09, -1.37558e+07, 2.46273e+09,
		6.06924e+07, -2.0196e+09, -2.42352e+09, -2.13135e+09, -1.19209e-07,
		-4.38199e+08, -3.22327e+08, -6.44655e+08, -1.67684e+09, 1.19209e-07,
		-3.42042e+08, -1.52654e+07, 2.48323e+08, -4.42499e+07, 2.39392e+08,
		-5.72912e+07, 2.32991e+08, -6.94093e+07, 2.26491e+08, -8.06834e+07,
		2.20134e+08, -3.58104e+08, -3.69365e+08, -8.74928e+08, -7.2193e+08,
		-3.79832e+08, -7.76661e+08, -2.94763e+10, -1.27133e+09, 2.87755e+09,
		-3.2463e+08, -1.06482e+09, 1.23533e+09, 4.36962e+07, 2.34988e+08,
		-1.33356e+07, 2.11186e+08, -7.16381e+07, -4.10647e+07, 4.42295e+09,
		-1.83533e+07, 2.46829e+09, 9.4849e+06, 1.82359e+09, -3.46495e+08,
		-2.28598e+09, -2.28598e+09, -2.00879e+09, 0, -4.12828e+08, -1.57449e+09,
		0, -3.21099e+08, -2.65115e+07, 2.32682e+08, -3.89477e+07, 2.27948e+08,
		-5.05889e+07, 2.22463e+08, -6.1467e+07, 2.16796e+08, -3.36026e+08,
		-9.02019e+08, 1.17772e+09, 2.33666e+09, -3.56226e+08, -3.65355e+08,
		-2.30131e+10, -1.28487e+09, 2.81986e+09, 1.18829e+09, 1.18829e+09,
		1.10913e+09, 1.10913e+09, 4.50548e+09, -4.16439e+07, 2.04407e+08,
		-4.03532e+07, -2.44589e+07, 2.48664e+09, -3.90035e+06, 1.81458e+09,
		2.06811e+07, -3.08994e+08, -2.05927e+09, 0, -1.807e+09, -1.4456e+09,
		-3.71077e+08, -8.43963e+08, -8.43963e+08, -2.86748e+08, -1.04509e+07,
		2.12779e+08, -2.09098e+07, 2.11773e+08, -3.09005e+07, 2.08492e+08,
		-2.49306e+08, -2.53287e+08, 1.07164e+09, 1.07164e+09, 2.13355e+09,
		-3.17525e+08, -3.25531e+08, -1.8604e+10, 1.32213e+09, 2.80212e+09,
		1.17234e+09, 1.17234e+09, 1.03909e+09, 1.03909e+09, 2.4983e+09,
		-2.64843e+07, 2.00145e+08, -2.78028e+07, -4.16478e+07, 4.54713e+09,
		-8.52002e+06, 1.8141e+09, 1.47999e+09, -2.93417e+08, -1.96459e+09,
		3.27432e+08, -1.72282e+09, -2.06739e+09, -3.53671e+08, -5.34735e+08,
		-8.02102e+08, -2.72474e+08, -9.35566e+06, 2.03423e+08, -1.87673e+07,
		2.0288e+08, -2.32991e+08, -2.37263e+08, -2.41005e+08, 1.08566e+09,
		1.08566e+09, 1.02243e+09, 1.02243e+09, -3.01452e+08, -3.08992e+08,
		-1.55817e+10, 1.31559e+09, 2.78959e+09, 1.15985e+09, 1.15985e+09,
		9.8765e+08, 9.8765e+08, 1.81554e+09, 1.81554e+09, 1.94826e+08,
		-1.69462e+07, -4.15299e+07, 4.58878e+09, -2.80446e+07, 2.5111e+09,
		1.474e+09, -2.79498e+08, -1.87926e+09, 6.2642e+08, -1.64739e+09,
		-1.97687e+09, -3.38077e+08, -2.54895e+08, -5.09789e+08, -2.59716e+08,
		-8.42838e+06, 1.94995e+08, -2.1771e+08, -2.22463e+08, -2.26491e+08,
		-2.3002e+08, 1.06698e+09, 1.06698e+09, 1.01325e+09, 1.01325e+09,
		1.9707e+09, -2.87089e+08, -2.94215e+08, -1.26074e+10, 1.31101e+09,
		2.78116e+09, 1.15007e+09, 1.15007e+09, 9.48653e+08, 9.48653e+08,
		1.47015e+09, 1.87359e+08, -7.63613e+06, -4.13236e+07, 4.63031e+09,
		2.52478e+09, 2.52478e+09, 1.81848e+09, 1.81848e+09, -2.6698e+08,
		-1.49468e+09, -1.19209e-07, -1.89524e+09, -1.89524e+09, -3.24021e+08,
		-2.48241e+08, -2.0288e+08, -2.08492e+08, -2.12986e+08, -2.16796e+08,
		-2.20134e+08, 1.05199e+09, 1.05199e+09, 9.92437e+08, 9.92437e+08,
		9.59847e+08, 9.59847e+08, -2.74174e+08, -2.80928e+08, -1.72076e+10,
		6.924e+07, 7.75842e+09, -3.73599e+09, 9.38305e+08, 5.83935e+08,
		2.93104e+09, -5.47012e+08, -3.35642e+08, -2.62495e+08, 9.15219e+07,
		4.13777e+09, -4.40012e+09, -1.45424e+09, 8.35125e+07, 305852, -7511.29,
		1.56112e+09, -1.33496e+09, -5.62443e+08, 6.3123e+07, 1.02718e+09,
		4.67532e+08, -3.46781e+08, -163879, -112625, -86484, -70549.6, -59784.1,
		-3.79832e+08, -3.56226e+08, -3.17525e+08, -3.01452e+08, -2.87089e+08,
		-2.74174e+08, -1.75368e+10, 1.11818e+09, 1.13852e+10, -7.01864e+09,
		1.80622e+09, 1.10343e+09, 3.24952e+09, -6.91634e+08, -3.44172e+08,
		-2.68913e+08, -8.52203e+09, -2.16157e+09, 8.13559e+07, 1.26108e+09,
		160139, -7265.33, -1.53021e+09, -6.37682e+08, 6.15457e+07, 6.29369e+08,
		5.0187e+08, 5.89e+07, -7087.95, 313364, -114968, -88238.6, -71951.7,
		-60951.6, -4.1353e+08, -7.76661e+08, -3.65355e+08, -3.25531e+08,
		-3.08992e+08, -2.94215e+08, -2.80928e+08, -9.33618e+09, 4.17875e+07,
		6.00778e+09, -2.54517e+09, 7.71891e+08, 8.40042e+08, 2.34997e+09,
		-4.62127e+08, 4.43418e+07, 3.12225e+09, -3.30096e+09, -1.10759e+09,
		3.95702e+07, 281564, -4480.31, 1.28523e+09, -1.0979e+09, 0, 4.60149e+08,
		1.16071e+09, 3.66307e+08, 3.66307e+08, 43629.5, 138657, 53800.8,
		92218.9, 215971, -9.93792e+09, 9.30321e+08, 8.61545e+09, -4.89533e+09,
		1.50069e+09, 1.26006e+09, 2.61944e+09, -5.85621e+08, -6.35844e+09,
		-1.63308e+09, 3.79796e+07, 9.18842e+08, 147469, -4300.6, -1.26462e+09,
		-5.31639e+08, 3.6353e+07, 9.20299e+08, 3.94774e+08, 3.94774e+08,
		-4120.88, 286044, 45387.9, 94326.6, 55508, 70266.1, -1.33205e+10,
		-1.50618e+09, 3.05997e+09, 1.72362e+08, -3.15598e+08, -2.98663e+08,
		-3.43584e+08, -6.59495e+08, -6.94176e+08, 3.61724e+08, 4.04099e+08,
		3.71825e+08, -7.40979e+07, -1.22055e+09, -2.58126e+09, -2.46662e+09,
		-2.1523e+09, -1.64895e+09, -345131, -266725, -180059, -88813.9, 5502.88,
		-36125.9, 176119, -87739.3, 14225.1, -8.65306e+09, -2.5101e+08,
		-1.84641e+09, -1.25315e+09, 4.59517e+09, 2.62203e+08, -4.23607e+06,
		2.89612e+08, 4.32587e+08, 4.26499e+08, -2.71324e+07, -5.32334e+07,
		6.80077e+08, -4.52863e+09, -2.09154e+09, -1.13184e+09, -2.75663e+08,
		2.88091e+08, -8.58961e+07, -3.71716e+08, -1.53439e+08, -217488, -247270,
		-221838, -170113, -98054.4, -68658.7, 69138.8, -3.3503e+07,
		-3.88443e+07, -4.73888e+07, -6.09942e+07, -1.84085e+10, -1.02364e+09,
		-3.33968e+09, 2.00104e+09, -1.36969e+09, 5.86142e+07, -6.06571e+08,
		7.17968e+07, -4.33121e+08, -2.86031e+08, -4.45644e+09, -3.44404e+09,
		-2.32497e+09, -1.14679e+09, -24635, -166686, -111180, -54269.4,
		-46271.7, -199907, -191029, -4.29788e+08, -8.76598e+08, -1.04047e+09,
		-9.46369e+08, -6.45861e+08, -6.8862e+08, -6.46271e+08, -7.61995e+08,
		-3.43212e+08, -3.27432e+08, -3.75852e+08, -2.98936e+08, -57616.5,
		-25026.4, -1.94726e+10, -1.94421e+09, -3.44823e+09, 7.97417e+08,
		-1.13623e+09, -5.271e+08, -3.80322e+08, -3.64244e+08, -2.80827e+09,
		-3.19284e+09, -2.86444e+09, -2.19656e+09, -163130, -80171.1, -350721,
		-161980, -87655.5, -21348.8, 59471.2, 26259.3, -8.06219e+08,
		-7.08071e+08, -6.33451e+08, -5.74676e+08, -4.8774e+08, -4.54596e+08,
		-4.2627e+08, -4.01759e+08, -5.05961e+08, -5.51303e+08, -5.27165e+08,
		-6.31747e+08, -12190.3, -1.1508e+10, -8.75483e+08, -2.22036e+09,
		1.39903e+09, -8.48729e+08, 3.32491e+07, -4.08225e+08, 5.9945e+07,
		-3.26758e+09, -2.54293e+09, -1.73007e+09, -8.66817e+08, -22689.3,
		-149202, -124650, -108645, -77857.9, -178334, -170770, -3.63128e+08,
		-7.19779e+08, -8.45681e+08, -7.64462e+08, -2.97211e+08, -2.75971e+08,
		-1.93396e+08, -52132.2, -22847.3, -1.36408e+10, -1.60273e+09,
		-2.22193e+09, 5.26639e+08, -9.29186e+08, -4.20021e+08, -2.37634e+08,
		-2.0198e+09, -2.29506e+09, -2.06153e+09, -1.58231e+09, -144743,
		-70502.1, -306587, -142816, -76944.2, -18007.4, 53382.5, 23206.3,
		-6.56377e+08, -5.73059e+08, -5.09837e+08, -4.60149e+08, -3.86903e+08,
		-3.59084e+08, -3.35368e+08, -3.14899e+08, -2.25057e+08, -2.13894e+08,
		-1.52937e+08, -11073.7};
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
					hPartials[startingIdx + index], 0.05, 1.0e-6);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::V] = 0;
	comp[Spec::I] = 2;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::I] = 0;
	comp[Spec::PerfectV] = 15;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 15);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 20);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::PerfectV] = 0;
	comp[Spec::FaultedI] = 12;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 10);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 15);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::FaultedI] = 0;
	comp[Spec::PerfectI] = 5;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 5);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.321142, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 7.2480, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 62637676, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 20736000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.5, 0.01);
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 6);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);

	comp[Spec::PerfectI] = 0;
	comp[Spec::FaultedV] = 20;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
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
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectV], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectV], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedV], 20);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedV], 21);
	BOOST_REQUIRE_EQUAL(lo[Spec::FaultedI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::FaultedI], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::PerfectI], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::PerfectI], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 1);
}

BOOST_AUTO_TEST_SUITE_END()
