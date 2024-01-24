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
	knownDFill[0] = {0, 9, 10, 11, 17, 18, 19, 20, 21, 22, 1, 2, 3, 12};
	knownDFill[1] = {1, 0, 9, 10, 17, 18, 19, 20, 21, 22, 2, 3, 12, 11};
	knownDFill[2] = {2, 0, 10, 1, 9, 17, 18, 19, 20, 21, 22, 3, 11, 12};
	knownDFill[3] = {3, 0, 11, 1, 10, 2, 9, 17, 18, 19, 20, 21, 22, 12};
	knownDFill[4] = {4, 9, 10, 11, 12, 17, 18, 19, 20, 21, 22, 5, 6, 7, 8};
	knownDFill[5] = {5, 4, 9, 10, 11, 17, 18, 19, 20, 21, 22, 6, 7, 8, 12};
	knownDFill[6] = {6, 4, 10, 5, 9, 17, 18, 19, 20, 21, 22, 7, 8, 12, 11};
	knownDFill[7] = {7, 4, 11, 5, 10, 6, 9, 17, 18, 19, 20, 21, 22, 8, 12};
	knownDFill[8] = {8, 4, 12, 5, 11, 6, 10, 7, 9, 17, 18, 19, 20, 21, 22};
	knownDFill[9] = {9, 0, 20, 1, 21, 2, 22, 4, 19, 5, 6, 7, 10, 11, 12, 13, 14,
		15, 16, 17, 18};
	knownDFill[10] = {10, 0, 19, 1, 20, 2, 21, 3, 22, 4, 18, 5, 6, 7, 8, 9, 11,
		12, 13, 14, 15, 16, 17};
	knownDFill[11] = {11, 0, 18, 1, 19, 2, 20, 3, 21, 4, 17, 5, 6, 7, 8, 9, 10,
		12, 13, 14, 15, 16, 22};
	knownDFill[12] = {12, 0, 17, 1, 18, 2, 19, 3, 20, 4, 5, 6, 7, 8, 9, 11, 10,
		13, 14, 15, 16, 21, 22};
	knownDFill[13] = {13, 9, 14, 10, 11, 12, 17, 22, 18, 21, 19, 20};
	knownDFill[14] = {14, 9, 10, 11, 12, 13, 17, 18, 22, 19, 21, 20};
	knownDFill[15] = {15, 9, 16, 10, 11, 12, 17, 22, 18, 21, 19, 20};
	knownDFill[16] = {16, 9, 10, 11, 12, 15, 17, 18, 22, 19, 21, 20};
	knownDFill[17] = {17, 0, 22, 1, 2, 3, 4, 21, 5, 6, 7, 8, 9, 18, 10, 19, 11,
		20, 12, 13, 15, 14, 16};
	knownDFill[18] = {
		18, 0, 1, 2, 3, 4, 22, 5, 6, 7, 8, 9, 19, 10, 20, 11, 21, 12, 17};
	knownDFill[19] = {19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 10, 21, 11, 22, 12,
		13, 15, 17, 18};
	knownDFill[20] = {20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 21, 10, 22, 11, 13, 15,
		12, 14, 16, 17, 19, 18};
	knownDFill[21] = {21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 22, 10, 13, 15, 11, 14,
		16, 12, 17, 20, 18, 19};
	knownDFill[22] = {22, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15, 10, 14, 16, 11,
		12, 17, 21, 18, 20, 19};

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
		network.getLeftSideRate(dConcs, 0, gridId), 7043583050, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-4.04847e+09, -1.64971e+06, 2.75808e+09,
		8.10829e+09, -1.39154e+10, -4.64112e+09, -6.56636e+08, 2.69417e+09,
		1.02258e+10, -7.68572e+10, -3.54494e+10, -2.29667e+10, -1.27326e+10,
		-4.48549e+09, -8.7327e+09, -3.46879e+09, -6.64667e+09, -1.22541e+09,
		-2.12621e+09, -1.96182e+09, -5.90354e+08, -1.24477e+08, 1.20908e+09, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-7.0437e+09, -1.9521e+09, -7.46054e+08,
		-1.36392e+08, 8712.16, 8360.78, 8062.23, -80436.6, -65709.3, -55747.7,
		297817, 160139, 112625, 1.78091e+09, -5.8623e+09, 3.73301e+09,
		-1.12494e+08, -8.93509e+08, 8035.56, 7746.79, -107467, -82617.9,
		-67456.9, -57206.2, 305852, 163879, 1.12251e+09, 2.0093e+09,
		-3.95005e+09, 1.9598e+09, 1.9598e+09, 3.8455e+09, -1.03758e+08, 7511.29,
		-160139, -110136, -84619.7, -69058.9, -58542, 313364, 9.02511e+08,
		9.02511e+08, -806684, 1.35014e+09, 1.35014e+09, 2.01602e+09,
		2.01602e+09, 3.94926e+09, 3.94926e+09, -313364, -163879, -112625,
		-86484, -70549.6, -59784.1, 1.58735e+09, -2.26534e+10, -2.42016e+09,
		1.13673e+08, -3.28098e+08, -1.91491e+09, -1.25051e+08, -2.22549e+08,
		-2.99122e+08, -3.61138e+08, -8.15975e+08, -8.31706e+08, 5.86249e+08,
		5.27139e+08, 4.77458e+08, 4.36815e+08, -1.32913e+10, 4.26622e+09,
		1.82809e+09, -1.37738e+09, -8.2978e+08, -6.711e+08, -1.55139e+08,
		-2.1271e+08, -6.55896e+08, -6.70397e+08, -6.83058e+08, 5.01398e+08,
		4.61942e+08, 4.25993e+08, 1.78091e+09, -1.05535e+10, 2.69838e+09,
		1.29133e+09, 4.21903e+09, 9.63549e+06, -6.15285e+07, -6.4187e+08,
		-5.45212e+08, -5.59593e+08, -5.71726e+08, -5.82322e+08, 4.3987e+08,
		4.12408e+08, 1.12251e+09, 2.0093e+09, -7.57613e+09, 2.17416e+09,
		3.07667e+09, 2.59113e+09, 2.59113e+09, 4.2094e+09, -1.06058e+07,
		-4.67617e+07, -4.61942e+08, -9.54917e+08, -4.89812e+08, -5.0024e+08,
		-5.09349e+08, 3.93108e+08, 9.02511e+08, -3.00503e+09, 5.3732e+09,
		3.50226e+09, 2.04353e+09, 2.04353e+09, 2.52956e+09, 2.52956e+09,
		4.22e+09, 4.22e+09, -3.93108e+08, -4.12408e+08, -4.25993e+08,
		-8.7363e+08, -4.45954e+08, -4.5394e+08, -9.25469e+10, -3.73293e+09,
		-2.19354e+09, -3.84544e+09, -2.34558e+09, -3.9492e+09, -2.46729e+09,
		-3.48964e+09, -1.21891e+09, -3.56314e+09, -3.63767e+09, -3.71065e+09,
		-2.66788e+09, -4.56704e+09, -2.78814e+09, -3.94926e+09, -4.04625e+09,
		-3.0297e+09, -3.07791e+09, -1.19014e+09, -1.69712e+09, -5.35083e+10,
		-1.9597e+09, -8.17784e+08, -2.01594e+09, -9.3708e+08, 69058.9,
		-1.03037e+09, 59784.1, -1.10633e+09, -1.94869e+09, 9.20152e+07,
		-1.95243e+09, -1.96997e+09, 5.0024e+08, 4.5394e+08, 6.09466e+09,
		-1.50422e+09, -1.49655e+09, -2.06776e+09, -2.11605e+09, -1.58681e+09,
		-1.61091e+09, -4.11028e+08, -3.45099e+10, -1.34999e+09, -3.53904e+08,
		107467, -4.59425e+08, 84619.7, -5.41617e+08, 70549.6, -6.08099e+08,
		-1.46286e+09, 5.00643e+08, -1.42645e+09, 5.45212e+08, 4.89812e+08,
		4.45954e+08, -1.53597e+09, 1.5747e+09, -1.09239e+09, -1.42212e+09,
		-1.45424e+09, -1.09152e+09, -1.10759e+09, -1.07389e+09, -2.44398e+10,
		289105, -1.2609e+08, 156132, -2.21187e+08, 110136, -2.96109e+08, 86484,
		-3.56846e+08, -1.91491e+09, 5.86249e+08, 5.27139e+08, 4.77458e+08,
		4.36815e+08, 1.30296e+08, 1.88709e+09, 3.37907e+09, -1.09264e+09,
		-1.11671e+09, -8.38734e+08, -8.50787e+08, -8.10701e+08, -8.25523e+08,
		-8.53208e+09, 9.69882e+07, 4.04625e+09, -2.06776e+09, -1.42212e+09,
		-1.09264e+09, -144483, 161369, 99573.2, 99573.2, 83874.5, 83874.5,
		-8.73326e+09, -4.04625e+09, -2.11605e+09, -1.45424e+09, -1.11671e+09,
		305852, 305852, 96669.4, 96669.4, 78367.2, 78367.2, 147519,
		-6.54705e+09, 4.82093e+07, 3.07791e+09, -1.58681e+09, -1.09152e+09,
		-8.38734e+08, -120194, 161369, 99573.2, 99573.2, 83874.5, 83874.5,
		-6.6472e+09, -3.07791e+09, -1.61091e+09, -1.10759e+09, -8.50787e+08,
		281564, 281564, 96669.4, 96669.4, 78367.2, 78367.2, 147519,
		-8.71822e+09, -233357, 6.82797e+08, -297817, -305852, -313364,
		1.04674e+08, 1.62637e+09, -4.8944e+08, -5.01398e+08, -4.3987e+08,
		-3.93108e+08, 1.43806e+08, 2.69309e+09, 1.22277e+08, 1.45667e+09,
		1.08228e+08, 1.03124e+09, 9.80725e+07, -303720, -280155, 2175.51,
		1348.72, -9.1049e+09, -151778, -156132, -160139, -163879, 8.20175e+07,
		1.65704e+09, -6.1708e+08, -1.05428e+09, -4.61942e+08, -4.12408e+08,
		1.00916e+08, 2.76937e+09, 9.06005e+07, 1.49655e+09, 8.29874e+07,
		1.05392e+09, 7.70413e+07, 566405, -9.34809e+09, -104562, -107467,
		-110136, -112625, -7.7658e+08, -6.38702e+08, -5.45212e+08, -9.54917e+08,
		-4.25993e+08, 8.03497e+07, 2.84941e+09, 7.40857e+07, 1.53052e+09,
		6.91443e+07, 1.07389e+09, 1.15769e+09, 1.09264e+09, 8.38734e+08,
		28817.5, 146344, -9.55011e+09, -80436.6, -82617.9, -84619.7, -86484,
		-7.97954e+08, -6.55896e+08, -5.59593e+08, -4.89812e+08, -8.7363e+08,
		6.78568e+07, 2.91739e+09, 6.35758e+07, 1.56033e+09, 1.48212e+09,
		1.42212e+09, 1.09152e+09, 1.17375e+09, 1.11671e+09, 8.50787e+08, 12109,
		175376, 314032, -9.31814e+09, -65709.3, -67456.9, -69058.9, -70549.6,
		-8.15975e+08, -6.70397e+08, -5.71726e+08, -5.0024e+08, -4.45954e+08,
		5.9317e+07, 2.9767e+09, 2.12389e+09, 2.06776e+09, 1.58681e+09,
		1.50771e+09, 1.45424e+09, 1.10759e+09, -8.10701e+08, 5671.81, 331016,
		26452.3, 68864.3, -9.49757e+09, -55747.7, -57206.2, -58542, -59784.1,
		-8.31706e+08, -6.83058e+08, -5.82322e+08, -5.09349e+08, -4.5394e+08,
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
	BOOST_REQUIRE_EQUAL(bounds.size(), 39);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 6);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 6);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 1);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 2, 39, 3, 40, 5, 6, 41, 7, 42, 9, 10, 11, 12, 14, 43,
		15, 44, 18, 45, 19, 46, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
		33, 34};
	knownDFill[1] = {1, 0, 5, 6, 41, 7, 42, 9, 10, 11, 12, 21, 22, 23, 24, 25,
		26, 31, 32, 33, 34, 2, 39, 27, 28, 29, 30};
	knownDFill[2] = {2, 0, 39, 27, 28, 29, 30, 1, 5, 31, 32, 33, 34, 6, 41, 9,
		10, 11, 12, 21, 22, 23, 24, 25, 26, 3, 40, 4};
	knownDFill[3] = {3, 0, 2, 39, 40, 1, 6, 41, 5, 9, 10, 11, 12, 31, 32, 33,
		34, 21, 22, 23, 24, 25, 26, 4, 27, 28, 29, 30};
	knownDFill[4] = {4, 0, 3, 40, 1, 7, 42, 2, 39, 6, 41, 31, 32, 33, 34, 5, 9,
		10, 11, 12, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
	knownDFill[5] = {5, 0, 9, 1, 2, 39, 3, 40, 6, 41, 7, 42, 10, 11, 12, 14, 43,
		15, 44, 18, 45, 19, 46, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
		33, 34};
	knownDFill[6] = {6, 0, 41, 31, 32, 33, 34, 1, 2, 39, 5, 9, 10, 11, 12, 15,
		44, 16, 19, 46, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 7, 42, 8};
	knownDFill[7] = {7, 0, 6, 41, 42, 1, 5, 9, 10, 11, 12, 31, 32, 33, 34, 16,
		20, 21, 22, 23, 24, 25, 26, 8};
	knownDFill[8] = {8, 0, 7, 42, 5, 6, 41, 31, 32, 33, 34, 9, 10, 11, 12, 21,
		22, 23, 24, 25, 26};
	knownDFill[9] = {9, 0, 23, 1, 24, 2, 39, 3, 40, 5, 6, 41, 7, 42, 10, 11, 12,
		13, 14, 43, 15, 44, 16, 17, 18, 45, 19, 46, 20, 21, 22, 25, 26, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
	knownDFill[10] = {10, 0, 22, 1, 23, 2, 39, 3, 40, 5, 6, 41, 7, 42, 9, 11,
		12, 13, 14, 43, 15, 44, 16, 17, 18, 45, 19, 46, 20, 21, 24, 25, 26, 27,
		28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
	knownDFill[11] = {11, 0, 21, 1, 22, 2, 39, 3, 40, 5, 6, 41, 7, 42, 9, 10,
		12, 13, 14, 43, 15, 44, 16, 17, 18, 45, 19, 46, 20, 23, 24, 25, 26, 27,
		28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
	knownDFill[12] = {12, 0, 1, 21, 2, 39, 26, 3, 40, 5, 6, 41, 7, 42, 9, 11,
		10, 13, 14, 43, 15, 44, 16, 17, 18, 45, 19, 46, 20, 22, 23, 24, 25, 27,
		28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
	knownDFill[13] = {13, 0, 14, 43, 5, 6, 41, 15, 44, 16, 9, 35, 10, 36, 11,
		12, 21, 22, 23, 24, 25, 26, 31, 33, 34};
	knownDFill[14] = {14, 0, 43, 15, 44, 5, 9, 10, 11, 12, 13, 23, 24, 25, 26,
		21, 22, 31, 32, 33, 34, 36, 35};
	knownDFill[15] = {15, 0, 44, 5, 6, 41, 9, 16, 10, 11, 12, 14, 43, 21, 22,
		23, 24, 25, 26, 31, 32, 33, 34, 36};
	knownDFill[16] = {16, 6, 41, 7, 42, 9, 10, 11, 12, 14, 43, 26, 15, 44, 21,
		22, 23, 24, 25};
	knownDFill[17] = {17, 0, 18, 45, 5, 6, 41, 19, 46, 20, 9, 37, 10, 38, 11,
		12, 21, 22, 23, 24, 25, 26, 31, 33, 34};
	knownDFill[18] = {18, 0, 45, 19, 46, 5, 9, 10, 11, 12, 17, 23, 24, 25, 26,
		21, 22, 31, 32, 33, 34, 38, 37};
	knownDFill[19] = {19, 0, 46, 5, 6, 41, 9, 20, 10, 11, 12, 18, 45, 21, 22,
		23, 24, 25, 26, 31, 32, 33, 34, 38};
	knownDFill[20] = {20, 6, 41, 7, 42, 9, 10, 11, 12, 18, 45, 26, 19, 46, 21,
		22, 23, 24, 25};
	knownDFill[21] = {21, 0, 25, 1, 26, 2, 39, 3, 40, 4, 5, 6, 41, 15, 44, 19,
		46, 7, 42, 16, 20, 8, 9, 22, 10, 23, 11, 24, 12, 13, 14, 43, 34, 17, 18,
		45, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38};
	knownDFill[22] = {22, 0, 26, 1, 2, 39, 3, 40, 4, 5, 6, 41, 15, 44, 19, 46,
		7, 42, 16, 20, 8, 9, 23, 10, 24, 11, 25, 12, 13, 14, 43, 33, 34, 17, 18,
		45, 21, 27, 28, 29, 30, 31, 32, 35, 36, 37, 38};
	knownDFill[23] = {23, 0, 1, 2, 39, 3, 40, 4, 5, 6, 41, 15, 44, 19, 46, 7,
		42, 16, 20, 8, 9, 24, 10, 25, 11, 26, 12, 13, 17, 14, 43, 32, 33, 34,
		18, 45, 21, 22, 27, 28, 29, 30, 31, 35, 36, 37, 38};
	knownDFill[24] = {24, 0, 1, 2, 39, 3, 40, 4, 5, 6, 41, 15, 44, 19, 46, 7,
		42, 16, 20, 8, 9, 25, 10, 26, 11, 13, 17, 12, 35, 37, 14, 43, 31, 32,
		33, 34, 18, 45, 21, 23, 22, 27, 28, 29, 30, 36, 38};
	knownDFill[25] = {25, 0, 1, 2, 39, 3, 40, 4, 5, 14, 43, 18, 45, 6, 41, 15,
		44, 19, 46, 7, 42, 8, 9, 26, 10, 13, 17, 11, 35, 37, 12, 36, 38, 31, 32,
		33, 21, 24, 22, 23, 27, 28, 29, 30, 34};
	knownDFill[26] = {26, 0, 14, 43, 18, 45, 1, 2, 39, 3, 40, 4, 5, 6, 41, 15,
		44, 16, 19, 46, 20, 7, 42, 8, 9, 13, 17, 10, 35, 37, 11, 36, 38, 12, 31,
		32, 34, 21, 25, 22, 24, 23, 27, 28, 29, 30, 33};
	knownDFill[27] = {27, 0, 1, 9, 2, 39, 24, 25, 26, 5, 6, 41, 10, 12, 11, 21,
		28, 22, 29, 23, 30, 31, 32, 33, 34};
	knownDFill[28] = {28, 0, 1, 10, 2, 39, 23, 24, 25, 26, 5, 6, 41, 9, 27, 11,
		12, 21, 29, 22, 30, 31, 32, 33, 34};
	knownDFill[29] = {29, 0, 1, 11, 2, 39, 22, 23, 24, 25, 26, 5, 6, 41, 9, 28,
		10, 27, 12, 21, 30, 31, 32, 33, 34};
	knownDFill[30] = {30, 0, 1, 12, 2, 39, 21, 22, 23, 24, 25, 3, 40, 26, 5, 6,
		41, 9, 29, 10, 28, 11, 27, 31, 32, 33, 34};
	knownDFill[31] = {31, 0, 10, 1, 2, 39, 5, 9, 6, 41, 24, 25, 26, 12, 11, 14,
		43, 15, 44, 18, 45, 19, 46, 21, 32, 22, 33, 23, 34, 27, 28, 29, 30};
	knownDFill[32] = {32, 0, 11, 1, 2, 39, 5, 10, 6, 41, 23, 24, 25, 26, 9, 31,
		12, 14, 43, 15, 44, 18, 45, 19, 46, 21, 33, 22, 34, 27, 28, 29, 30};
	knownDFill[33] = {33, 0, 12, 1, 2, 39, 5, 11, 6, 41, 22, 23, 24, 25, 26, 9,
		32, 10, 31, 14, 43, 15, 44, 18, 45, 19, 46, 21, 34, 27, 28, 29, 30};
	knownDFill[34] = {34, 0, 5, 1, 2, 39, 12, 6, 41, 21, 22, 23, 24, 25, 7, 42,
		26, 9, 33, 10, 32, 11, 31, 14, 43, 15, 44, 18, 45, 19, 46, 27, 28, 29,
		30};
	knownDFill[35] = {35, 0, 14, 43, 5, 6, 41, 15, 44, 16, 9, 36, 10, 11, 12,
		13, 21, 32, 33, 34, 22, 26, 23, 25, 24};
	knownDFill[36] = {36, 0, 14, 43, 6, 41, 16, 9, 10, 11, 12, 13, 22, 15, 44,
		31, 32, 33, 34, 21, 35, 23, 26, 24, 25};
	knownDFill[37] = {37, 0, 18, 45, 5, 6, 41, 19, 46, 20, 9, 38, 10, 11, 12,
		17, 21, 32, 33, 34, 22, 26, 23, 25, 24};
	knownDFill[38] = {38, 0, 18, 45, 6, 41, 20, 9, 10, 11, 12, 17, 22, 19, 46,
		31, 32, 33, 34, 21, 37, 23, 26, 24, 25};
	knownDFill[39] = {39, 0, 2, 27, 28, 29, 30, 1, 5, 31, 32, 33, 34, 6, 41, 9,
		10, 11, 12, 21, 22, 23, 24, 25, 26, 3, 40, 4};
	knownDFill[40] = {40, 0, 2, 39, 3, 1, 6, 41, 5, 9, 10, 11, 12, 31, 32, 33,
		34, 21, 22, 23, 24, 25, 26, 4, 27, 28, 29, 30};
	knownDFill[41] = {41, 0, 6, 31, 32, 33, 34, 1, 2, 39, 5, 9, 10, 11, 12, 15,
		44, 16, 19, 46, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 7, 42, 8};
	knownDFill[42] = {42, 0, 6, 41, 7, 1, 5, 9, 10, 11, 12, 31, 32, 33, 34, 16,
		20, 21, 22, 23, 24, 25, 26, 8};
	knownDFill[43] = {43, 0, 14, 15, 44, 5, 9, 10, 11, 12, 13, 23, 24, 25, 26,
		21, 22, 31, 32, 33, 34, 36, 35};
	knownDFill[44] = {44, 0, 15, 5, 6, 41, 9, 16, 10, 11, 12, 14, 43, 21, 22,
		23, 24, 25, 26, 31, 32, 33, 34, 36};
	knownDFill[45] = {45, 0, 18, 19, 46, 5, 9, 10, 11, 12, 17, 23, 24, 25, 26,
		21, 22, 31, 32, 33, 34, 38, 37};
	knownDFill[46] = {46, 0, 19, 5, 6, 41, 9, 20, 10, 11, 12, 18, 45, 21, 22,
		23, 24, 25, 26, 31, 32, 33, 34, 38};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 1452);
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
		network.getLeftSideRate(dConcs, 0, gridId), 63192846675, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-4.40818e+10, -1.09636e+10, -5.15591e+09,
		8.06425e+09, 2.71298e+10, -4.05186e+10, -1.09852e+10, 9.0959e+09,
		2.59987e+10, -2.79011e+11, -1.32572e+11, -8.68898e+10, -6.32018e+10,
		3.32575e+09, 1.97022e+09, -2.47657e+09, -1.25068e+10, 2.79577e+09,
		1.60301e+09, -2.29447e+09, -8.96595e+09, -5.38964e+09, -6.80409e+09,
		-7.13308e+09, -5.18967e+09, -5.49413e+09, 1.91072e+09, -8.82862e+09,
		-6.93208e+09, -4.86288e+09, -3.26817e+09, -3.59807e+10, -2.60566e+10,
		-1.64166e+10, -1.19992e+10, 4.12316e+09, -5.10793e+09, 3.54652e+09,
		-3.5819e+09, -1.17491e+10, -7.08746e+09, -2.15596e+10, -9.92027e+09,
		-1.21177e+10, -1.65828e+10, -9.02504e+09, -1.24088e+10, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-5.28373e+10, -9.33235e+08,
		-5.33884e+09, 0, -2.28264e+09, 3.42396e+09, -1.01014e+09, -6.59906e+09,
		-7.71512e+08, -2.67814e+09, 4.0172e+09, -2.42016e+09, 1.13673e+08,
		-3.28098e+08, -1.91491e+09, -4.27107e+09, 2.13554e+09, -4.56528e+09,
		2.28264e+09, -3.60102e+09, 1.80051e+09, -3.76617e+09, 1.88309e+09,
		-1.25051e+08, -2.22549e+08, -2.99122e+08, -3.61138e+08, -4.12715e+08,
		-1.21746e+09, -9.56868e+08, -9.78472e+08, -9.98528e+08, -1.01736e+09,
		-9.84007e+08, -9.75587e+08, -9.7506e+08, -9.78688e+08, -1.39589e+10,
		-9.33235e+08, -7.75244e+08, -2.22143e+09, 0, -3.12017e+08, 6.24033e+08,
		-1.9521e+09, -7.46054e+08, -1.36392e+08, 7.42286e+08, 8712.16, 8360.78,
		8062.23, 7802.01, -141424, -55747.7, -6.67579e+08, -5.89157e+08,
		-5.2671e+08, -4.75189e+08, 139797, -215512, 297817, 160139, 112625,
		88238.6, -7.1904e+09, -4.91077e+08, -4.22673e+09, 9.29126e+08,
		1.12363e+09, 1.44492e+09, 2.15666e+09, 6.06776e+08, -7.26897e+07,
		-2.05415e+08, 9.49109e+07, 2.13613e+08, 2.0351e+08, 1.03219e+08,
		-2.06438e+08, -1.803e+09, -1.42854e+09, -1.00686e+09, -5.57166e+08,
		-5590.15, -2778.61, 0, 2758.5, 5502.88, 32532.7, 444131, -352238,
		14225.1, -2.45361e+09, 1.28132e+09, 6.78053e+09, 3.69282e+09,
		-1.52763e+09, 4.44286e+08, 1.39181e+09, -9.60257e+08, 1.0713e+09,
		-2.16545e+08, 1.29868e+09, 1.54204e+09, 1.44027e+09, 4.47805e+08,
		2.52915e+08, 2.17921e+08, 2.96887e+08, 150004, 77810, 27092.8, -19226.1,
		-64380.8, -68658.7, 180801, 2.8263e+08, 2.10968e+08, 2.42533e+08,
		3.91539e+08, -975132, 0, 1.22658e+10, 7.64003e+09, -3.12017e+08,
		3.12017e+08, -6.24033e+08, 2.04726e+09, 2.67127e+09, 1.40302e+09,
		1.93156e+09, 2.27255e+09, 1.33283e+09, 5.97059e+08, 0, -9.41595e+08,
		1.42355e+10, 4.92692e+09, 1.68497e+09, 0, -382527, -198266, -135464,
		-103565, -84184, -71125.3, 1.38062e+09, 9.4965e+08, 4.8822e+08, 0,
		-5.2136e+10, 2.66982e+09, 1.82809e+09, -7.75244e+08, -4.4149e+09, 0,
		-9.41595e+08, 1.88319e+09, -5.43811e+09, -1.14447e+09, -1.16602e+09,
		2.33204e+09, -1.37738e+09, -8.2978e+08, 1.01967e+07, -3.53192e+09,
		1.76596e+09, -4.70797e+09, 0, -2.95297e+09, 1.47648e+09, -3.85745e+09,
		0, -6.711e+08, -1.55139e+08, -2.1271e+08, -2.60766e+08, -1.04975e+09,
		-6.83058e+08, -7.94246e+08, -8.11585e+08, -8.27656e+08, -8.42722e+08,
		-8.80665e+08, -8.55218e+08, -8.44257e+08, -8.40196e+08, -1.71827e+10,
		-5.24753e+08, -2.5513e+09, 7.77875e+08, 1.54628e+09, 1.9559e+09,
		2.74298e+09, -4.44286e+08, 1.03219e+08, -2.06438e+08, 3.90608e+08,
		-1.81388e+09, -1.54981e+09, -1.14023e+09, -5.93144e+08, -2.39864e+09,
		6.54174e+08, -2.24268e+08, -1.895e+09, 5.16819e+08, -1.7574e+08,
		1.92403e+07, 2.02937e+07, 0, -4.32342e+07, -1.10556e+08, 6.1285e+07,
		-4.60207e+08, -1.8993e+08, 0, 9.91734e+07, 1.00198e+09, -3.98369e+08,
		4.22371e+07, -4.29611e+09, 1.67636e+09, 8.79378e+09, 3.52868e+09,
		-4.09118e+08, 6.24033e+07, 1.47385e+09, -1.04359e+08, 1.51351e+09,
		1.89518e+09, 1.87466e+09, 6.82815e+08, 4.37009e+08, 6.51523e+08,
		8.46484e+08, -1.60272e+08, -1.23487e+08, 8.56752e+07, 8.91832e+07,
		4.0029e+07, -6.47016e+07, -2.27197e+08, -2.73044e+08, 1.98951e+08,
		-1.20594e+09, 0, 1.3048e+10, 6.62636e+09, -1.16602e+09, 2.23394e+09,
		6.27223e+09, 3.4126e+09, 2.13007e+09, 1.01176e+09, 0, 1.38633e+10,
		5.11866e+09, 1.85845e+09, 0, -1.87359e+08, -1.94826e+08, -2.00145e+08,
		-2.04407e+08, -2.0802e+08, -2.11186e+08, -2.947e+11, -3.48964e+09,
		-1.21891e+09, -3.73293e+09, -2.19354e+09, -2.19213e+10, 0, -2.37259e+10,
		0, -3.56314e+09, -2.22358e+10, 0, -2.31054e+10, 0, -2.66788e+09,
		-4.56704e+09, -2.78814e+09, -3.94926e+09, -2.19213e+10, 0, -2.37259e+10,
		0, -4.93931e+09, -3.0297e+09, -1.61809e+10, 0, -1.69656e+10, 0,
		-3.47348e+09, -1.19014e+09, -1.69712e+09, -2.34558e+09, -2.46729e+09,
		-3.84544e+09, -3.9492e+09, -4.04625e+09, -4.13777e+09, -3.63767e+09,
		-3.71065e+09, -4.24241e+09, -4.27206e+09, -4.04625e+09, -4.13777e+09,
		-3.07791e+09, -3.12225e+09, -1.50631e+11, -1.94869e+09, 9.20152e+07,
		-1.9597e+09, -8.17784e+08, -1.14203e+10, 0, -9.85384e+09, 4.92692e+09,
		-1.95243e+09, -1.26862e+10, 0, -1.02373e+10, 5.11866e+09, 6.09466e+09,
		-1.50422e+09, -1.49655e+09, -2.06776e+09, -1.14203e+10, 0, -1.23173e+10,
		0, -2.56008e+09, -1.58681e+09, -8.45026e+09, 0, -8.8426e+09, 0,
		-1.8087e+09, -4.11028e+08, -9.3708e+08, -1.03037e+09, -1.10633e+09,
		-2.01594e+09, -2.06769e+09, -2.11599e+09, -2.16157e+09, -1.96997e+09,
		-1.99343e+09, -2.01959e+09, -2.46363e+09, -2.11605e+09, -2.16157e+09,
		-1.61091e+09, -1.63308e+09, -9.88434e+10, -1.46286e+09, 5.00643e+08,
		-1.34999e+09, -3.53904e+08, -7.82931e+09, 0, -5.0549e+09, 5.0549e+09,
		-1.42645e+09, -9.44531e+09, 0, -5.57536e+09, 5.57536e+09, -1.53597e+09,
		1.5747e+09, -1.09239e+09, -1.42212e+09, -7.82931e+09, 0, -8.42484e+09,
		0, -1.74915e+09, -1.09152e+09, -5.80176e+09, 0, -6.06332e+09, 0,
		-1.23945e+09, -4.59425e+08, -5.41617e+08, -6.08099e+08, -6.63462e+08,
		-1.38754e+09, -1.42203e+09, -1.45417e+09, -1.48445e+09, -1.41719e+09,
		-1.41961e+09, -1.42795e+09, -1.43955e+09, -1.45424e+09, -1.48451e+09,
		-1.10759e+09, -1.12237e+09, -7.49265e+10, -1.91491e+09, -1.03833e+09,
		-1.2609e+08, -6.00104e+09, -128165, -1.21134e+09, -2.57854e+09,
		3.8678e+09, -1.18446e+09, -7.42533e+09, -7.71512e+08, -3.00807e+09,
		4.5121e+09, 1.30296e+08, 1.88709e+09, 3.37907e+09, -1.09264e+09,
		-6.00111e+09, 0, -6.44634e+09, 0, -1.33727e+09, -8.38734e+08,
		-4.45178e+09, 2.38419e-07, -4.64795e+09, 0, -9.49679e+08, -2.21187e+08,
		-2.96109e+08, -3.56846e+08, -4.07369e+08, -1.06663e+09, -1.09253e+09,
		-1.11662e+09, -1.13929e+09, -1.15141e+09, -1.13866e+09, -1.13543e+09,
		-1.13742e+09, -1.11671e+09, -1.13936e+09, -8.50787e+08, -8.61872e+08,
		-8.53257e+09, 0, 5.47435e+09, -4.6422e+09, 8.8298e+08, 3.83154e+09,
		1.61991e+09, 4.479e+09, 3.39002e+08, 0, 9.69882e+07, 4.04625e+09,
		9.38122e+07, 2.16157e+09, -2.98798e+09, -1.09264e+09, -144483, -60565.5,
		-26261.9, -745.235, 30514.3, 102827, 1.51503e+09, -6.33451e+08, 0,
		-7.30188e+09, 2.94205e+07, 6.34904e+09, 7.15466e+09, -4.23196e+09,
		5.88402e+08, -7.21818e+07, -3.58783e+07, 0, 3.56187e+07, 64471.4,
		-77976.8, -46887.7, 2692.36, 5096.4, -139634, -110634, 6.64457e+08,
		5.83119e+08, 3.94627e+08, -2.16945e+08, 152603, 98664.2, -1.08378e+10,
		-4.56528e+08, 6.24327e+09, -9.41595e+08, -2.39864e+09, 6.54174e+08,
		1.9369e+09, 2.11716e+09, 1.00471e+09, 3.4983e+08, -2.48254e+08, 412768,
		326260, -16770.4, 100577, 119424, 111542, 91958.2, 37823.4,
		-8.06219e+08, -7.08071e+08, -6.33451e+08, -2.2987e+08, 12190.3,
		-1.25085e+10, -2.24268e+09, 1.12134e+09, -1.60272e+09, 8.01362e+08,
		-4.93931e+09, -2.56008e+09, -1.74915e+09, -1.33727e+09, 64082.7, 128165,
		192248, 869834, 663484, 1.10248e+06, 381567, 130493, 0, -81217.6,
		-6.54749e+09, 0, 4.31773e+09, -3.48279e+09, 7.38242e+08, 3.02345e+09,
		1.28338e+09, 3.55408e+09, 2.42906e+08, 0, 4.82093e+07, 3.07791e+09,
		4.62756e+07, 1.63308e+09, -2.25188e+09, -8.38734e+08, -120194, -47896,
		-17565.7, 5927.09, 35953.8, 107434, 1.25705e+09, -5.09837e+08, 0,
		-5.72918e+09, 1.6515e+07, 4.85907e+09, 5.4549e+09, -2.98212e+09,
		4.76194e+08, -3.13868e+07, -1.56934e+07, 0, 1.56934e+07, 59388.4,
		-67668, -39085.6, 7332.55, 10103.9, -122419, -96768.4, 5.36243e+08,
		4.68036e+08, 3.1432e+08, -1.76559e+08, 138336, 90181.5, -8.38457e+09,
		-3.76617e+08, 4.59518e+09, -7.71491e+08, -1.895e+09, 5.16819e+08,
		1.37332e+09, 1.49426e+09, 7.15443e+08, 2.4789e+08, -1.819e+08, 367039,
		289487, -8750.75, 91321.6, 106865, 99293.9, 81532.2, 33927.1,
		-6.56377e+08, -5.73059e+08, -5.09837e+08, -1.8406e+08, 11073.7,
		-8.96744e+09, -1.7574e+09, 8.78701e+08, -1.23487e+09, 6.17437e+08,
		-3.47348e+09, -1.8087e+09, -1.23945e+09, -9.49679e+08, 57133.4, 114267,
		171400, 749158, 567893, 946015, 328713, 112698, 0, -70374.3,
		-1.21651e+10, 1.04674e+08, 1.62637e+09, -233357, 6.82797e+08,
		-1.69771e+06, 0, -1.83746e+06, 0, -382527, -4.8944e+08, -2.66192e+09,
		1.95165e+09, -5.46982e+08, 1.09029e+09, -4.32259e+08, 8.61365e+08,
		-4.95685e+08, -1.4188e+09, -4.00681e+08, -3.08719e+08, -1.87359e+08,
		1.43806e+08, 2.69309e+09, 1.22277e+08, 1.45667e+09, 1.08228e+08,
		1.03124e+09, 9.80725e+07, -303720, 5.40668e+08, -1.08473e+09,
		-1.34006e+09, -280155, 4.39896e+08, -8.82797e+08, -297817, -305852,
		-313364, -320452, -5.01398e+08, -4.3987e+08, -3.93108e+08, -311188,
		-319124, -284695, -289335, -1.26969e+10, 8.20175e+07, 1.65704e+09,
		-151778, -884453, 0, -953918, 0, -198266, -6.1708e+08, -2.73862e+09,
		9.75827e+08, -5.45908e+08, 5.45526e+08, -4.3134e+08, 4.31011e+08,
		-5.46049e+08, -7.094e+08, 0, 0, -1.94826e+08, 1.00916e+08, 2.76937e+09,
		9.06005e+07, 1.49655e+09, 8.29874e+07, 1.05392e+09, 7.70413e+07,
		-160139, 1.13854e+09, -1.73648e+09, -1.49826e+09, -3.73434e+08, -147469,
		9.29408e+08, -1.41899e+09, 566405, -156132, -160139, -163879, -167404,
		-1.05428e+09, -4.61942e+08, -163879, -167404, -149709, -151770,
		-1.30453e+10, -7.7658e+08, -104562, -606344, 0, -652465, 0, -135464,
		-6.38702e+08, -1.81697e+09, -1.95165e+09, -391479, -1.0899e+09, -338095,
		-8.61027e+08, -5.81792e+08, 0, 4.00681e+08, 3.08719e+08, -2.00145e+08,
		8.03497e+07, 2.84941e+09, 7.40857e+07, 1.53052e+09, 6.91443e+07,
		1.07374e+09, 1.15769e+09, 1.09253e+09, 8.38633e+08, 1.80523e+09,
		-1.92989e+09, -2.17042e+09, -4.25993e+08, 5.98255e+08, 1.47874e+09,
		-1.58697e+09, 28817.5, 146344, -107467, -110136, -112625, -114968,
		-5.45212e+08, -112625, -114968, -102934, -104307, -1.33323e+10,
		-7.97954e+08, -80436.6, -464758, -2.91038e-11, -499240, 0, -103565,
		-6.55896e+08, 9.14396e+07, -5.85496e+09, 1.09009e+09, -3.27057e+09,
		8.61192e+08, -2.58384e+09, -6.10388e+08, 7.094e+08, 8.01362e+08,
		6.17437e+08, -2.04407e+08, 6.78568e+07, 2.91726e+09, 6.35758e+07,
		1.5602e+09, 1.48212e+09, 1.42203e+09, 1.09144e+09, 1.17375e+09,
		1.11662e+09, 8.50708e+08, 2.56289e+09, -1.63908e+09, -1.94563e+09,
		-4.89812e+08, 2.12223e+08, 1.57239e+09, 2.10739e+09, -1.36474e+09,
		12109, 175376, 314032, -82617.9, -84619.7, -86484, -88238.6, -88238.6,
		-80097.7, -1.35761e+10, -8.15975e+08, -65709.3, -378574, 0, -406088, 0,
		-84184, -2.29162e+09, 2.90359e+09, -1.92642e+09, 2.40431e+09,
		-1.61621e+09, 2.98235e+09, -9.75827e+09, 2.72564e+09, -5.45129e+09,
		2.15334e+09, -4.30668e+09, -1.34401e+09, 0, -2.0802e+08, 5.9317e+07,
		2.97659e+09, 2.12389e+09, 2.06769e+09, 1.58674e+09, 1.50771e+09,
		1.45417e+09, 1.10753e+09, 1.19054e+09, 1.13929e+09, 8.61807e+08,
		-5.71726e+08, 7.15262e+08, 1.72575e+09, 5671.81, 195410, 26452.3,
		68864.3, -67456.9, -69058.9, -70549.6, -71951.7, -8.0652e+08,
		-1.30299e+10, -2.79973e+09, 4.57458e+09, -4.75254e+09, 3.70618e+09,
		-3.77038e+09, -55747.7, -320414, -1.45519e-11, -343293, 0, -71125.3,
		-6.83058e+08, 4.54045e+09, -3.18018e+09, 3.84555e+09, -2.23964e+09,
		-5.60669e+08, 3.04424e+09, -1.78166e+09, -4.3935e+08, -1.36522e+09, 0,
		-2.11186e+08, 4.00227e+09, 3.94921e+09, 3.02965e+09, 2.16663e+09,
		2.11599e+09, 1.61086e+09, 1.53299e+09, 1.48445e+09, 1.12232e+09,
		-2.9161e+09, 8.03718e+08, 1.92165e+09, -1.44519e+09, 2605.09, 210516,
		15301.8, 80782.5, 223772, -57206.2, -58542, -59784.1, -60951.6,
		-4.5394e+08, -1.46893e+10, -9.56868e+08, 3.73301e+09, -1.12494e+08,
		232749, -261618, -175570, -67456.9, 6876.53, -7.94246e+08, -2.30104e+09,
		0, -8.93509e+08, 5.57253e+07, 6.21653e+08, 8035.56, 305852, 7746.79,
		163879, 7501.03, 114968, -6.83498e+08, -6.02874e+08, -5.41353e+08,
		-4.92737e+08, -1.36405e+10, -9.78472e+08, 1.9598e+09, -1.07957e+08,
		354018, -271407, -231405, -84619.7, 6655.83, 69623.5, -8.11585e+08,
		-1.8993e+09, 9.4965e+08, -1.03758e+08, 3.8455e+09, -5.19604e+08,
		-1.90125e+08, 7511.29, 313364, 7265.33, 167404, -6.98001e+08,
		-6.15354e+08, -5.52318e+08, -5.02524e+08, -1.29726e+10, -9.98528e+08,
		1.35014e+09, -1.04102e+08, 530909, -271170, -340769, -112625, 6467.66,
		80879.9, 132464, -8.27656e+08, -1.46466e+09, 1.46466e+09, -9.69882e+07,
		3.94926e+09, -1.00029e+08, 2.01602e+09, 4.70639e+08, 7087.95, 320452,
		-7.11424e+08, -6.26891e+08, -5.62443e+08, -5.11553e+08, -1.27418e+10,
		-8.41224e+07, 1.97186e+09, -1.00742e+08, 806367, -611592, -659993,
		-167404, 6300.64, 97664.7, 155192, 68658.7, -137317, -129610,
		-8.42722e+08, -9.91734e+08, 1.4876e+09, -9.15219e+07, 4.04625e+09,
		-9.38122e+07, 2.06776e+09, -9.68557e+07, 1.38765e+09, -7.23993e+08,
		-6.37682e+08, -5.71905e+08, -5.19983e+08, -4.88752e+10, 1.18723e+09,
		1.29133e+09, -6.67579e+08, -3.78758e+09, 0, 2.83697e+09, 9.63549e+06,
		-4.5505e+09, -1.12344e+09, -9.31636e+08, -5.71726e+08, -1.96566e+08,
		-5.5603e+08, 4.6904e+07, -3.03006e+09, 1.51503e+09, -4.0311e+09, 0,
		-2.5141e+09, 1.25705e+09, -3.28188e+09, 0, -6.15285e+07, -7.89173e+08,
		-6.4187e+08, -7.65315e+08, -1.59706e+08, -7.53506e+08, -6.83498e+08,
		-6.98001e+08, -7.11424e+08, -7.23993e+08, -3.98245e+10, 7.2111e+08,
		1.16725e+09, -5.89157e+08, -2.66566e+09, 1.33283e+09, 1.27397e+09,
		9.7461e+07, -2.75959e+09, 1.41699e+09, -1.31831e+09, -4.89812e+08,
		-1.20883e+08, 2.62162e+08, -1.06058e+07, 2.98035e+09, -7.13606e+08,
		-2.66566e+09, 1.33283e+09, -3.54036e+09, 0, -2.19635e+09, 1.09817e+09,
		-2.86529e+09, 0, -4.67617e+07, -7.20298e+08, -8.85077e+07, -6.98203e+08,
		-6.02874e+08, -6.15354e+08, -6.26891e+08, -6.37682e+08, -3.40728e+10,
		3.96132e+09, 1.93002e+09, -5.2671e+08, -1.79118e+09, 1.79118e+09,
		7.73281e+08, 1.69627e+08, -1.18219e+09, 3.11765e+09, -7.64966e+08,
		-4.25993e+08, -5.01587e+08, 3.1276e+08, 7.03327e+08, -2.24086e+07,
		3.1066e+09, 5.60325e+07, 1.35184e+09, -2.38824e+09, 1.19412e+09,
		-3.16725e+09, 0, -1.95518e+09, 9.77589e+08, -2.54919e+09, 0,
		-3.68048e+07, -6.66122e+08, -5.41353e+08, -5.52318e+08, -5.62443e+08,
		-5.71905e+08, -2.76708e+10, 2.14445e+08, 2.13178e+09, -4.75189e+08,
		-1.08473e+09, 1.62709e+09, 2.30028e+08, -1.38141e+08, 3.01783e+09,
		-6.93522e+08, -3.73434e+08, -2.21125e+07, 3.48956e+08, 3.31552e+08,
		2.73044e+08, -5.46088e+08, -6.8341e+08, -2.96523e+07, 3.21999e+09,
		3.00402e+07, 1.42204e+09, 1.12486e+08, 8.23387e+08, -2.16945e+09,
		1.08473e+09, -2.2987e+09, 1.14935e+09, -1.76559e+09, 8.82797e+08,
		-1.8406e+09, 9.20299e+08, -4.92737e+08, -5.02524e+08, -5.11553e+08,
		-5.19983e+08, -8.73406e+09, 1.06777e+09, 7.0009e+09, -5.25102e+09,
		1.76596e+09, 3.28639e+09, 3.81601e+09, 4.64192e+09, 1.76642e+09,
		5.60669e+08, 9.15219e+07, 4.13777e+09, -4.40012e+09, -1.45424e+09,
		8.35125e+07, 305852, -7511.29, -7.08071e+08, 0, 5.74676e+08, -67209.1,
		36885.4, -34257.4, 7817.6, 61035.4, -8.92404e+09, 2.13554e+09,
		1.05022e+10, -8.78461e+09, 5.60669e+08, 5.60669e+08, 1.12134e+09,
		-8.52203e+09, -2.16157e+09, 8.13559e+07, 1.26108e+09, 160139, -7265.33,
		2.72242e+09, -1.74583e+09, -8.06219e+08, 0, 6.33451e+08, 1.14935e+09,
		-7087.95, 313364, -40120.7, 13896, -20435.7, -4148.82, -6.64794e+09,
		9.00255e+08, 5.37926e+09, -3.80222e+09, 1.47648e+09, 2.59276e+09,
		3.01478e+09, 3.69646e+09, 1.35882e+09, 4.3935e+08, 4.43418e+07,
		3.12225e+09, -3.30096e+09, -1.10759e+09, 3.95702e+07, 281564, -4480.31,
		-5.73059e+08, 0, 4.60149e+08, -53039.9, 41987.6, -24566.4, 13851.8,
		68451.9, -6.74033e+09, 1.80051e+09, 7.87721e+09, -6.37182e+09,
		4.3935e+08, 4.3935e+08, 8.78701e+08, -6.35844e+09, -1.63308e+09,
		3.79796e+07, 9.18842e+08, 147469, -4300.6, 2.19942e+09, -1.42566e+09,
		-6.56377e+08, 0, 5.09837e+08, 9.20299e+08, -4120.88, 286044, -29459.6,
		19479, -12294.8, 2463.24, -1.29149e+10, -1.46881e+09, 2.85714e+09,
		-3.15598e+08, -2.98663e+08, -3.43584e+08, -6.59495e+08, -7.40979e+07,
		-8.66144e+08, -6.12189e+08, -8.3208e+07, 4.15298e+07, -1.34349e+07,
		-1.03219e+08, 2.06438e+08, -2.58126e+09, -2.46662e+09, -2.1523e+09,
		-1.64895e+09, -345131, -266725, -180059, -88813.9, 5502.88, -36125.9,
		176119, -87739.3, 14225.1, -8.65306e+09, -2.42974e+08, -1.84641e+09,
		-1.25315e+09, 4.59517e+09, 4.44286e+08, 4.80128e+08, 1.87251e+08,
		6.94661e+08, -4.52863e+09, -2.09154e+09, -1.13184e+09, -2.75663e+08,
		3.09711e+08, -5.58634e+07, -2.26414e+08, -2.05274e+08, -217488, -247270,
		-221838, -170113, -98054.4, -68658.7, 69138.8, 1.77577e+08, -6.5751e+07,
		-2.30445e+08, -2.55742e+08, -2.30165e+10, -1.98569e+09, 3.62939e+09,
		-1.55523e+09, -5.24577e+08, -5.21879e+08, -7.55719e+08, -4.44286e+08,
		-1.03219e+08, 2.06438e+08, -1.59506e+09, -2.63327e+09, -2.77086e+09,
		-2.64724e+09, -2.21791e+09, -2.39864e+09, 6.54174e+08, -2.24268e+08,
		-1.895e+09, 5.16819e+08, -1.7574e+08, -3.17978e+08, -5.18691e+08,
		-5.66122e+08, -4.36894e+08, -1.10556e+08, -3.66061e+08, -4.60207e+08,
		-1.8993e+08, -9.7644e+07, -9.91734e+07, 1.99185e+08, 2.60935e+08,
		4.22371e+07, -1.15062e+10, -2.56235e+08, -1.76434e+09, -1.00309e+09,
		4.9242e+09, -6.24033e+07, 1.00744e+09, -4.51673e+09, -2.28794e+09,
		-1.31316e+09, -2.8917e+08, 4.54717e+08, -1.16838e+08, -5.37213e+08,
		-3.32167e+08, -1.60272e+08, -1.23487e+08, -1.55342e+08, -3.32043e+08,
		-4.39005e+08, -4.42773e+08, -3.10405e+08, -2.73044e+08, -5.0902e+06,
		-1.17248e+10, 2.57684e+08, -1.36888e+09, 2.11598e+09, -1.13982e+09,
		5.88402e+08, -4.45644e+09, -3.44404e+09, -2.32497e+09, -1.14679e+09,
		-24635, -166686, -123966, -84623.1, -64798.9, -199907, -191029,
		-1.41763e+08, -5.49795e+08, -6.82239e+08, -2.16945e+08, -57616.5,
		-25026.4, -1.59119e+10, -4.56528e+08, -1.78082e+09, -9.41595e+08,
		-2.39864e+09, 6.54174e+08, -2.80827e+09, 1.11014e+09, -3.19284e+09,
		-2.86444e+09, -2.19656e+09, -163130, -80171.1, -350721, -161980,
		-87655.5, -21348.8, 59471.2, 26259.3, -8.06219e+08, -7.08071e+08,
		-6.33451e+08, -2.2987e+08, -12190.3, -8.92029e+09, 2.04824e+08,
		-9.31014e+08, 1.49106e+09, -6.64669e+08, 4.76194e+08, -3.26758e+09,
		-2.54293e+09, -1.73007e+09, -8.66817e+08, -22689.3, -149202, -111865,
		-78291.6, -59330.7, -178334, -170770, -1.20134e+08, -4.48858e+08,
		-5.52404e+08, -1.76559e+08, -52132.2, -22847.3, -1.19693e+10,
		-3.76617e+08, -1.21988e+09, -7.71491e+08, -1.895e+09, 5.16819e+08,
		-2.0198e+09, 7.80597e+08, -2.29506e+09, -2.06153e+09, -1.58231e+09,
		-144743, -70502.1, -306587, -142816, -76944.2, -18007.4, 53382.5,
		23206.3, -6.56377e+08, -5.73059e+08, -5.09837e+08, -1.8406e+08,
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
