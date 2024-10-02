#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/AlloyReactionNetwork.h>
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
 * This suite is responsible for testing the Alloy network.
 */
BOOST_AUTO_TEST_SUITE(AlloyNetwork_testSuite)

BOOST_AUTO_TEST_CASE(fullyRefined)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
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

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 22);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 22);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.011664, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 22);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 6);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 6);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 1);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 8, 9, 10, 16, 17, 18, 19, 20, 21, 1, 2, 3, 11};
	knownDFill[1] = {1, 0, 8, 9, 16, 17, 18, 19, 20, 21, 2, 3, 11, 10};
	knownDFill[2] = {2, 0, 9, 1, 8, 16, 17, 18, 19, 20, 21, 3, 10, 11};
	knownDFill[3] = {3, 0, 10, 1, 9, 2, 8, 16, 17, 18, 19, 20, 21, 11};
	knownDFill[4] = {4, 8, 9, 10, 16, 17, 18, 19, 20, 21, 5, 6, 7, 11};
	knownDFill[5] = {5, 4, 8, 9, 16, 17, 18, 19, 20, 21, 6, 7, 11, 10};
	knownDFill[6] = {6, 4, 9, 5, 8, 16, 17, 18, 19, 20, 21, 7, 10, 11};
	knownDFill[7] = {7, 4, 10, 5, 9, 6, 8, 16, 17, 18, 19, 20, 21, 11};
	knownDFill[8] = {
		8, 0, 19, 1, 20, 2, 21, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
	knownDFill[9] = {9, 0, 18, 1, 19, 2, 20, 3, 21, 4, 5, 6, 7, 8, 10, 11, 12,
		13, 14, 15, 16, 17};
	knownDFill[10] = {10, 0, 17, 1, 18, 2, 19, 3, 20, 4, 5, 6, 7, 8, 9, 11, 12,
		13, 14, 15, 16, 21};
	knownDFill[11] = {11, 0, 16, 1, 17, 2, 18, 3, 19, 4, 5, 6, 7, 8, 9, 10, 12,
		13, 14, 15, 20, 21};
	knownDFill[12] = {12, 8, 13, 9, 10, 11, 16, 21, 17, 20, 18, 19};
	knownDFill[13] = {13, 8, 9, 10, 11, 12, 16, 17, 18, 21, 19, 20};
	knownDFill[14] = {14, 8, 15, 9, 10, 11, 16, 17, 18, 19, 21, 20};
	knownDFill[15] = {15, 8, 9, 10, 11, 14, 16, 17, 21, 18, 19, 20};
	knownDFill[16] = {16, 0, 21, 1, 2, 3, 4, 5, 6, 7, 8, 17, 9, 18, 10, 19, 11,
		20, 12, 14, 13, 15};
	knownDFill[17] = {
		17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 18, 9, 19, 10, 20, 11, 21, 16};
	knownDFill[18] = {
		18, 0, 1, 2, 3, 4, 5, 6, 7, 8, 19, 9, 20, 10, 21, 11, 12, 14, 16, 17};
	knownDFill[19] = {19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 20, 9, 21, 10, 12, 14, 11,
		13, 15, 16, 17, 18};
	knownDFill[20] = {20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 9, 12, 13, 15, 10, 11,
		14, 16, 17, 19, 18};
	knownDFill[21] = {21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 14, 9, 13, 15, 10, 11,
		16, 17, 20, 18, 19};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 372);
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
		network.getTotalAtomConcentration(dConcs, Spec::Perfect, 3), 26.0,
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
		8.10829e+09, -8.32109e+09, -2.82787e+09, 9.97476e+08, 7.01855e+09,
		-7.15215e+10, -3.06887e+10, -1.96578e+10, -1.08177e+10, -4.48549e+09,
		-8.73270e+09, -3.46879e+09, -6.64667e+09, -7.43837e+08, -1.68109e+09,
		-7.07786e+08, 6.44415e+08, 6.91498e+08, 2.04079e+09, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-7.04370e+09, -1.95210e+09,
		-7.46054e+08, -1.36392e+08, 8.71216e+03, 8.36078e+03, 8.06223e+03,
		-8.04366e+04, -6.57093e+04, -5.57477e+04, 2.97817e+05, 1.60139e+05,
		1.12625e+05, 1.78091e+09, -5.86230e+09, 3.73301e+09, -1.12494e+08,
		-8.93509e+08, 8.03556e+03, 7.74679e+03, -1.07467e+05, -8.26179e+04,
		-6.74569e+04, -5.72062e+04, 3.05852e+05, 1.63879e+05, 1.12251e+09,
		2.00930e+09, -3.95005e+09, 1.95980e+09, 1.95980e+09, 3.84550e+09,
		-1.03758e+08, 7.51129e+03, -1.60139e+05, -1.10136e+05, -8.46197e+04,
		-6.90589e+04, -5.85420e+04, 3.13364e+05, 9.02511e+08, 9.02511e+08,
		-8.06684e+05, 1.35014e+09, 1.35014e+09, 2.01602e+09, 2.01602e+09,
		3.94926e+09, 3.94926e+09, -3.13364e+05, -1.63879e+05, -1.12625e+05,
		-8.64840e+04, -7.05496e+04, -5.97841e+04, 1.58735e+09, -1.27051e+10,
		-2.43812e+09, -1.37738e+09, -8.29780e+08, -8.48509e+07, -1.55139e+08,
		-2.12710e+08, -6.55896e+08, -6.70397e+08, -6.83058e+08, 5.01398e+08,
		4.61942e+08, 4.25993e+08, 1.78091e+09, -1.00263e+10, 4.21903e+09,
		9.63549e+06, -1.40705e+09, -6.15285e+07, -1.14731e+08, -5.45212e+08,
		-5.59593e+08, -5.71726e+08, -5.82322e+08, 4.39870e+08, 4.12408e+08,
		1.12251e+09, 2.00930e+09, -7.09868e+09, 2.59113e+09, 2.59113e+09,
		4.20940e+09, -1.06058e+07, -4.67617e+07, -4.61942e+08, -4.77458e+08,
		-4.89812e+08, -5.00240e+08, -5.09349e+08, 3.93108e+08, 9.02511e+08,
		9.02511e+08, -2.56822e+09, 2.04353e+09, 2.04353e+09, 2.52956e+09,
		2.52956e+09, 4.22000e+09, 4.22000e+09, -3.93108e+08, -4.12408e+08,
		-4.25993e+08, -4.36815e+08, -4.45954e+08, -4.53940e+08, 1.58735e+09,
		-8.64347e+10, -3.73293e+09, -2.19354e+09, -3.84544e+09, -2.34558e+09,
		-3.94920e+09, -2.46729e+09, -3.56314e+09, -3.63767e+09, -3.71065e+09,
		-2.66788e+09, -2.72098e+09, -2.78814e+09, -3.94926e+09, -4.04625e+09,
		-3.02970e+09, -3.07791e+09, -1.19014e+09, -1.69712e+09, -1.99549e+09,
		-4.51858e+10, -1.95970e+09, -8.17784e+08, -2.01594e+09, -9.37080e+08,
		6.90589e+04, -1.03037e+09, 5.97841e+04, -1.10633e+09, -1.95243e+09,
		-1.96997e+09, 5.00240e+08, 4.53940e+08, 6.09466e+09, -1.50422e+09,
		-1.49655e+09, -2.06776e+09, -2.11605e+09, -1.58681e+09, -1.61091e+09,
		-4.11028e+08, -6.57673e+08, -3.04897e+10, -1.34999e+09, -3.53904e+08,
		1.07467e+05, -4.59425e+08, 8.46197e+04, -5.41617e+08, 7.05496e+04,
		-6.08099e+08, -1.42645e+09, 5.45212e+08, 4.89812e+08, 4.45954e+08,
		3.10083e+08, 1.57470e+09, -1.09239e+09, -1.42212e+09, -1.45424e+09,
		-1.09152e+09, -1.10759e+09, -2.10657e+08, -1.07389e+09, -2.25249e+10,
		2.89105e+05, -1.26090e+08, 1.56132e+05, -2.21187e+08, 1.10136e+05,
		-2.96109e+08, 8.64840e+04, -3.56846e+08, 5.86249e+08, 5.27139e+08,
		4.77458e+08, 4.36815e+08, 1.30296e+08, 3.37907e+09, 1.88709e+09,
		-1.09264e+09, -1.11671e+09, -8.38734e+08, -8.50787e+08, -8.10701e+08,
		-8.25523e+08, -8.53208e+09, 9.69882e+07, 4.04625e+09, -2.06776e+09,
		-1.42212e+09, -1.09264e+09, -1.44483e+05, 1.61369e+05, 9.95732e+04,
		9.95732e+04, 8.38745e+04, 8.38745e+04, -8.73326e+09, -4.04625e+09,
		-2.11605e+09, -1.45424e+09, -1.11671e+09, 3.05852e+05, 3.05852e+05,
		9.66694e+04, 7.83672e+04, 9.66694e+04, 1.47519e+05, 7.83672e+04,
		-6.54705e+09, 4.82093e+07, 3.07791e+09, -1.58681e+09, -1.09152e+09,
		-8.38734e+08, -1.20194e+05, 9.95732e+04, 8.38745e+04, 8.38745e+04,
		1.61369e+05, 9.95732e+04, -6.64720e+09, -3.07791e+09, -1.61091e+09,
		-1.10759e+09, -8.50787e+08, 2.81564e+05, 2.81564e+05, 9.66694e+04,
		9.66694e+04, 7.83672e+04, 1.47519e+05, 7.83672e+04, -7.42067e+09,
		-2.33357e+05, 6.82797e+08, -2.97817e+05, -3.05852e+05, -3.13364e+05,
		9.68089e+07, -5.01398e+08, -4.39870e+08, -3.93108e+08, 1.43806e+08,
		2.69309e+09, 1.22277e+08, 1.45667e+09, 1.08228e+08, 1.03124e+09,
		9.80725e+07, 8.10391e+08, -3.03720e+05, -2.80155e+05, 2.17551e+03,
		1.34872e+03, -7.82808e+09, -1.51778e+05, -1.56132e+05, -1.60139e+05,
		-1.63879e+05, -6.17080e+08, -5.27139e+08, -4.61942e+08, -4.12408e+08,
		1.00916e+08, 2.76937e+09, 9.06005e+07, 1.49655e+09, 8.29874e+07,
		1.05392e+09, 7.70413e+07, 8.25330e+08, 5.66405e+05, -8.09405e+09,
		-1.04562e+05, -1.07467e+05, -1.10136e+05, -1.12625e+05, -6.38702e+08,
		-5.45212e+08, -4.77458e+08, -4.25993e+08, 8.03497e+07, 2.84941e+09,
		7.40857e+07, 1.53052e+09, 6.91443e+07, 1.07389e+09, 1.15769e+09,
		1.09264e+09, 8.38734e+08, 2.88175e+04, 1.46344e+05, -8.31534e+09,
		-8.04366e+04, -8.26179e+04, -8.46197e+04, -8.64840e+04, -6.55896e+08,
		-5.59593e+08, -4.89812e+08, -4.36815e+08, 6.78568e+07, 2.91739e+09,
		6.35758e+07, 1.56033e+09, 1.48212e+09, 1.42212e+09, 1.09152e+09,
		1.17375e+09, 1.11671e+09, 8.50787e+08, 1.21090e+04, 3.14032e+05,
		1.75376e+05, -8.50217e+09, -6.57093e+04, -6.74569e+04, -6.90589e+04,
		-7.05496e+04, -6.70397e+08, -5.71726e+08, -5.00240e+08, -4.45954e+08,
		5.93170e+07, 2.97670e+09, 2.12389e+09, 2.06776e+09, 1.45424e+09,
		1.10759e+09, 1.50771e+09, -8.10701e+08, 1.58681e+09, 5.67181e+03,
		2.64523e+04, 3.31016e+05, 6.88643e+04, -8.66586e+09, -5.57477e+04,
		-5.72062e+04, -5.85420e+04, -5.97841e+04, -6.83058e+08, -5.82322e+08,
		-5.09349e+08, -4.53940e+08, 4.00227e+09, 3.94926e+09, 3.02970e+09,
		2.16663e+09, 2.11605e+09, 1.61091e+09, -1.07389e+09, -8.25523e+08,
		2.60509e+03, 1.53018e+04, 3.25344e+05, 3.73467e+05, 2.08641e+05};

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
				BOOST_REQUIRE(iter != row.end());
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 16);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 9);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 14);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 4);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 12);
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
	xolotl::options::ConfOptions opts;
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

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 38);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 46);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.011664, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.36, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 38);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 6);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 6);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 1);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;

	knownDFill[0] = {0, 4, 5, 40, 6, 41, 8, 9, 10, 11, 20, 21, 22, 23, 24, 25,
		30, 31, 32, 33, 1, 38, 26, 27, 28, 29};
	knownDFill[1] = {1, 0, 4, 30, 31, 32, 33, 38, 5, 40, 8, 9, 10, 11, 20, 21,
		22, 23, 24, 25, 2, 39, 3, 26, 27, 28, 29};
	knownDFill[2] = {2, 0, 5, 40, 1, 38, 4, 8, 9, 10, 11, 30, 31, 32, 33, 39,
		20, 21, 22, 23, 24, 25, 3, 26, 27, 28, 29};
	knownDFill[3] = {3, 0, 6, 41, 1, 38, 5, 40, 30, 31, 32, 33, 2, 39, 4, 8, 9,
		10, 11, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
	knownDFill[4] = {4, 0, 1, 38, 2, 39, 5, 40, 6, 41, 8, 9, 10, 11, 13, 42, 14,
		43, 17, 44, 18, 45, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
		33};
	knownDFill[5] = {5, 0, 40, 1, 38, 4, 30, 31, 32, 33, 8, 9, 10, 11, 14, 43,
		15, 18, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 6, 41, 7};
	knownDFill[6] = {6, 0, 41, 4, 5, 40, 8, 9, 10, 11, 30, 31, 32, 33, 15, 19,
		20, 21, 22, 23, 24, 25, 7};
	knownDFill[7] = {7, 4, 6, 41, 5, 40, 30, 31, 32, 33, 8, 9, 10, 11, 20, 21,
		22, 23, 24, 25};
	knownDFill[8] = {8, 0, 23, 1, 38, 2, 39, 4, 5, 40, 6, 41, 9, 10, 11, 12, 13,
		42, 14, 43, 15, 16, 17, 44, 18, 45, 19, 20, 21, 22, 24, 25, 26, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37};
	knownDFill[9] = {9, 0, 22, 1, 38, 2, 39, 4, 5, 40, 6, 41, 8, 10, 11, 12, 13,
		42, 14, 43, 15, 16, 17, 44, 18, 45, 19, 20, 21, 23, 24, 25, 26, 27, 28,
		29, 30, 31, 32, 33, 34, 35, 36, 37};
	knownDFill[10] = {10, 0, 21, 1, 38, 2, 39, 4, 5, 40, 6, 41, 8, 9, 11, 12,
		13, 42, 14, 43, 15, 16, 17, 44, 18, 45, 19, 20, 22, 23, 24, 25, 26, 27,
		28, 29, 30, 31, 32, 33, 34, 35, 36, 37};
	knownDFill[11] = {11, 0, 20, 1, 38, 25, 2, 39, 4, 5, 40, 6, 41, 8, 10, 9,
		12, 13, 42, 14, 43, 15, 16, 17, 44, 18, 45, 19, 21, 22, 23, 24, 26, 27,
		28, 29, 30, 31, 32, 33, 34, 35, 36, 37};
	knownDFill[12] = {12, 4, 13, 42, 5, 40, 14, 43, 15, 8, 34, 9, 35, 10, 11,
		20, 21, 22, 23, 24, 25, 30, 32, 33};
	knownDFill[13] = {13, 4, 42, 14, 43, 8, 9, 10, 11, 12, 22, 23, 24, 25, 20,
		21, 30, 31, 32, 33, 35, 34};
	knownDFill[14] = {14, 4, 43, 5, 40, 8, 15, 9, 10, 11, 13, 42, 20, 21, 22,
		23, 24, 25, 30, 31, 32, 33, 35};
	knownDFill[15] = {
		15, 5, 40, 6, 41, 8, 9, 10, 11, 13, 42, 25, 14, 43, 20, 21, 22, 23, 24};
	knownDFill[16] = {16, 4, 17, 44, 5, 40, 18, 45, 19, 8, 36, 9, 37, 10, 11,
		20, 21, 22, 23, 24, 25, 30, 32, 33};
	knownDFill[17] = {17, 4, 44, 18, 45, 8, 9, 10, 11, 16, 22, 23, 24, 25, 20,
		21, 30, 31, 32, 33, 37, 36};
	knownDFill[18] = {18, 4, 45, 5, 40, 8, 19, 9, 10, 11, 17, 44, 20, 21, 22,
		23, 24, 25, 30, 31, 32, 33, 37};
	knownDFill[19] = {
		19, 5, 40, 6, 41, 8, 9, 10, 11, 17, 44, 25, 18, 45, 20, 21, 22, 23, 24};
	knownDFill[20] = {20, 0, 25, 1, 38, 2, 39, 3, 4, 5, 40, 14, 43, 18, 45, 6,
		41, 15, 19, 7, 8, 21, 9, 22, 10, 23, 11, 24, 12, 13, 42, 33, 16, 17, 44,
		26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37};
	knownDFill[21] = {21, 0, 1, 38, 2, 39, 3, 4, 5, 40, 14, 43, 18, 45, 6, 41,
		15, 19, 7, 8, 22, 9, 23, 10, 24, 11, 25, 12, 13, 42, 32, 33, 16, 17, 44,
		20, 26, 27, 28, 29, 30, 31, 34, 35, 36, 37};
	knownDFill[22] = {22, 0, 1, 38, 2, 39, 3, 4, 5, 40, 14, 43, 18, 45, 6, 41,
		15, 19, 7, 8, 23, 9, 24, 10, 25, 11, 12, 16, 13, 42, 31, 32, 33, 17, 44,
		20, 21, 26, 27, 28, 29, 30, 34, 35, 36, 37};
	knownDFill[23] = {23, 0, 1, 38, 2, 39, 3, 4, 5, 40, 14, 43, 18, 45, 6, 41,
		15, 19, 7, 8, 24, 9, 25, 10, 12, 16, 11, 34, 36, 13, 42, 30, 31, 32, 33,
		17, 44, 20, 22, 21, 26, 27, 28, 29, 35, 37};
	knownDFill[24] = {24, 0, 1, 38, 2, 39, 3, 4, 13, 42, 17, 44, 5, 40, 14, 43,
		18, 45, 6, 41, 7, 8, 25, 9, 12, 16, 10, 34, 36, 11, 35, 37, 30, 31, 32,
		20, 23, 21, 22, 26, 27, 28, 29, 33};
	knownDFill[25] = {25, 0, 1, 38, 2, 39, 3, 4, 13, 42, 17, 44, 5, 40, 14, 43,
		15, 18, 45, 19, 6, 41, 7, 8, 12, 16, 9, 34, 36, 10, 35, 37, 11, 30, 31,
		33, 20, 24, 21, 23, 22, 26, 27, 28, 29, 32};
	knownDFill[26] = {26, 0, 8, 1, 38, 23, 24, 25, 4, 5, 40, 9, 11, 10, 20, 27,
		21, 28, 22, 29, 30, 31, 32, 33};
	knownDFill[27] = {27, 0, 9, 1, 38, 22, 23, 24, 25, 4, 5, 40, 8, 26, 10, 11,
		20, 28, 21, 29, 30, 31, 32, 33};
	knownDFill[28] = {28, 0, 10, 1, 38, 21, 22, 23, 24, 25, 4, 5, 40, 8, 27, 9,
		26, 11, 20, 29, 30, 31, 32, 33};
	knownDFill[29] = {29, 0, 11, 1, 38, 20, 21, 22, 23, 24, 2, 39, 25, 4, 5, 40,
		8, 28, 9, 27, 10, 26, 30, 31, 32, 33};
	knownDFill[30] = {30, 0, 1, 38, 4, 8, 5, 40, 23, 24, 25, 9, 11, 10, 13, 42,
		14, 43, 17, 44, 18, 45, 20, 31, 21, 32, 22, 33, 26, 27, 28, 29};
	knownDFill[31] = {31, 0, 1, 38, 4, 9, 5, 40, 22, 23, 24, 25, 8, 30, 10, 11,
		13, 42, 14, 43, 17, 44, 18, 45, 20, 32, 21, 33, 26, 27, 28, 29};
	knownDFill[32] = {32, 0, 1, 38, 4, 10, 5, 40, 21, 22, 23, 24, 25, 8, 31, 9,
		30, 11, 13, 42, 14, 43, 17, 44, 18, 45, 20, 33, 26, 27, 28, 29};
	knownDFill[33] = {33, 0, 1, 38, 4, 11, 5, 40, 20, 21, 22, 23, 24, 6, 41, 25,
		8, 32, 9, 31, 10, 30, 13, 42, 14, 43, 17, 44, 18, 45, 26, 27, 28, 29};
	knownDFill[34] = {34, 4, 13, 42, 5, 40, 14, 43, 15, 8, 35, 9, 10, 11, 12,
		20, 31, 32, 33, 21, 25, 22, 24, 23};
	knownDFill[35] = {35, 5, 40, 15, 8, 13, 42, 9, 10, 11, 12, 21, 14, 43, 30,
		31, 32, 33, 20, 34, 22, 25, 23, 24};
	knownDFill[36] = {36, 4, 17, 44, 5, 40, 18, 45, 19, 8, 37, 9, 10, 11, 16,
		20, 31, 32, 33, 21, 25, 22, 24, 23};
	knownDFill[37] = {37, 5, 40, 19, 8, 17, 44, 9, 10, 11, 16, 21, 18, 45, 30,
		31, 32, 33, 20, 36, 22, 25, 23, 24};
	knownDFill[38] = {38, 0, 4, 30, 31, 32, 33, 1, 5, 40, 8, 9, 10, 11, 20, 21,
		22, 23, 24, 25, 2, 39, 3, 26, 27, 28, 29};
	knownDFill[39] = {39, 0, 5, 40, 1, 38, 4, 8, 9, 10, 11, 30, 31, 32, 33, 2,
		20, 21, 22, 23, 24, 25, 3, 26, 27, 28, 29};
	knownDFill[40] = {40, 0, 5, 1, 38, 4, 30, 31, 32, 33, 8, 9, 10, 11, 14, 43,
		15, 18, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 6, 41, 7};
	knownDFill[41] = {41, 0, 6, 4, 5, 40, 8, 9, 10, 11, 30, 31, 32, 33, 15, 19,
		20, 21, 22, 23, 24, 25, 7};
	knownDFill[42] = {42, 4, 13, 14, 43, 8, 9, 10, 11, 12, 22, 23, 24, 25, 20,
		21, 30, 31, 32, 33, 35, 34};
	knownDFill[43] = {43, 4, 14, 5, 40, 8, 15, 9, 10, 11, 13, 42, 20, 21, 22,
		23, 24, 25, 30, 31, 32, 33, 35};
	knownDFill[44] = {44, 4, 17, 18, 45, 8, 9, 10, 11, 16, 22, 23, 24, 25, 20,
		21, 30, 31, 32, 33, 37, 36};
	knownDFill[45] = {45, 4, 18, 5, 40, 8, 19, 9, 10, 11, 17, 44, 20, 21, 22,
		23, 24, 25, 30, 31, 32, 33, 37};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 1371);
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
		network.getTotalAtomConcentration(dConcs, Spec::Perfect, 3), 200.0,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::Frank, 6), 9.377118,
		0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::Faulted, 1), 0.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 13649652815, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	std::vector<double> knownFluxes = {-1.00304e+10, -4.66483e+09, 6.78293e+09,
		2.71298e+10, -4.26021e+10, -1.05376e+10, 7.41954e+09, 2.59987e+10,
		-2.73675e+11, -1.27811e+11, -8.35809e+10, -6.12869e+10, 3.32575e+09,
		1.94080e+09, -2.02004e+09, -1.25068e+10, 2.79577e+09, 1.58649e+09,
		-1.91785e+09, -8.96595e+09, -4.90807e+09, -6.35897e+09, -5.87904e+09,
		-3.95491e+09, -4.27489e+09, 4.32470e+09, -7.87175e+09, -5.95361e+09,
		-3.86435e+09, -3.18405e+09, -3.66407e+10, -2.63002e+10, -1.82120e+10,
		-1.18104e+10, 3.05539e+09, -7.24346e+09, 2.64626e+09, -5.38241e+09,
		-1.02803e+10, -6.84449e+09, -1.94967e+10, -9.66404e+09, -1.23754e+10,
		-1.61263e+10, -9.22986e+09, -1.20322e+10, 0};
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i], knownFluxes[i], 0.01);
	}

	// Check the partials computation
	std::vector<double> knownPartials = {-1.30256e+10, -7.75244e+08,
		-2.22143e+09, 0.00000e+00, -3.12017e+08, 6.24033e+08, -1.95210e+09,
		-7.46054e+08, -1.36392e+08, 7.42286e+08, 8.71216e+03, 8.36078e+03,
		8.06223e+03, 7.80201e+03, -1.41424e+05, -5.57477e+04, -6.67579e+08,
		-5.89157e+08, -5.26710e+08, -4.75189e+08, 1.39797e+05, -2.15512e+05,
		2.97817e+05, 1.60139e+05, 1.12625e+05, 8.82386e+04, -6.33619e+09,
		6.06776e+08, -7.26897e+07, -2.05415e+08, 9.49109e+07, 2.13613e+08,
		2.03510e+08, -3.79962e+09, 1.03219e+08, -2.06438e+08, -1.80300e+09,
		-1.42854e+09, -1.00686e+09, -5.57166e+08, -5.59015e+03, -2.77861e+03,
		0.00000e+00, 2.75850e+03, 5.50288e+03, 3.25327e+04, 4.44131e+05,
		-3.52238e+05, 1.42251e+04, 7.37752e+08, 9.27938e+08, 1.24522e+09,
		1.95319e+09, -2.22535e+09, 4.44286e+08, 1.39181e+09, -9.60257e+08,
		5.92631e+09, 3.26572e+09, 1.07130e+09, -2.16545e+08, 1.29868e+09,
		1.54204e+09, 1.44027e+09, 4.47805e+08, 2.52915e+08, 2.17921e+08,
		2.96887e+08, -1.75589e+09, 1.50004e+05, 7.78100e+04, 2.70928e+04,
		-1.92261e+04, -6.43808e+04, -6.86587e+04, 1.80801e+05, 2.82630e+08,
		2.10968e+08, 2.42533e+08, 3.91539e+08, -9.75132e+05, -3.12017e+08,
		3.12017e+08, -6.24033e+08, 2.04726e+09, 2.67127e+09, 1.40302e+09,
		1.93156e+09, 2.27255e+09, 1.33283e+09, 5.97059e+08, 0.00000e+00,
		1.11245e+10, 8.78135e+09, -9.41595e+08, 1.42355e+10, 4.92692e+09,
		1.68497e+09, 0.00000e+00, -3.82527e+05, -1.98266e+05, -1.35464e+05,
		-1.03565e+05, -8.41840e+04, -7.11253e+04, 1.38062e+09, 9.49650e+08,
		4.88220e+08, 0.00000e+00, -4.99534e+10, -7.75244e+08, -4.41490e+09,
		0.00000e+00, -9.41595e+08, 1.88319e+09, -5.43811e+09, -1.14447e+09,
		-1.16602e+09, 2.33204e+09, -2.43812e+09, -1.37738e+09, -8.29780e+08,
		1.01967e+07, -3.53192e+09, 1.76596e+09, -4.70797e+09, 0.00000e+00,
		-2.95297e+09, 1.47648e+09, -3.85745e+09, 0.00000e+00, -8.48509e+07,
		-1.55139e+08, -2.12710e+08, -2.60766e+08, -1.04975e+09, -6.83058e+08,
		-7.94246e+08, -8.11585e+08, -8.27656e+08, -8.42722e+08, -8.80665e+08,
		-8.55218e+08, -8.44257e+08, -8.40196e+08, -1.59880e+10, -4.44286e+08,
		-2.14682e+09, 1.03219e+08, -2.06438e+08, 3.90608e+08, 4.75646e+08,
		1.25567e+09, 1.67353e+09, 2.46660e+09, -1.81388e+09, -1.54981e+09,
		-1.14023e+09, -5.93144e+08, -2.39864e+09, 6.54174e+08, -2.24268e+08,
		-1.89500e+09, 5.16819e+08, -1.75740e+08, 1.92403e+07, 2.02937e+07,
		0.00000e+00, -4.32342e+07, -1.10556e+08, -1.58661e+07, -4.60207e+08,
		-1.89930e+08, 0.00000e+00, 9.91734e+07, 1.00198e+09, -3.98369e+08,
		4.22371e+07, -4.02830e+09, 6.24033e+07, -6.76932e+08, 1.47385e+09,
		7.67621e+09, 2.96990e+09, -1.04359e+08, 1.51351e+09, 1.89518e+09,
		1.87466e+09, 6.82815e+08, 4.37009e+08, 6.51523e+08, 8.46484e+08,
		-1.60272e+08, -1.23487e+08, 8.56752e+07, 8.91832e+07, 4.00290e+07,
		-6.47016e+07, -2.27197e+08, -2.73044e+08, 1.98951e+08, -1.20594e+09,
		-1.16602e+09, 1.17089e+10, 7.96542e+09, 2.23394e+09, 6.27223e+09,
		3.41260e+09, 2.13007e+09, 1.01176e+09, 0.00000e+00, 1.38633e+10,
		5.11866e+09, 1.85845e+09, 0.00000e+00, -1.87359e+08, -1.94826e+08,
		-2.00145e+08, -2.04407e+08, -2.08020e+08, -2.11186e+08, -2.88588e+11,
		-3.73293e+09, -2.19354e+09, -2.19213e+10, 0.00000e+00, -2.37259e+10,
		0.00000e+00, -3.56314e+09, -2.22358e+10, 0.00000e+00, -2.31054e+10,
		0.00000e+00, -2.66788e+09, -2.72098e+09, -2.78814e+09, -3.94926e+09,
		-2.19213e+10, 0.00000e+00, -2.37259e+10, 0.00000e+00, -4.93931e+09,
		-3.02970e+09, -1.61809e+10, 0.00000e+00, -1.69656e+10, 0.00000e+00,
		-3.47348e+09, -1.19014e+09, -1.69712e+09, -1.99549e+09, -2.34558e+09,
		-2.46729e+09, -3.84544e+09, -3.94920e+09, -4.04625e+09, -4.13777e+09,
		-3.63767e+09, -3.71065e+09, -4.24241e+09, -4.27206e+09, -4.04625e+09,
		-4.13777e+09, -3.07791e+09, -3.12225e+09, -1.42308e+11, -1.95970e+09,
		-8.17784e+08, -1.14203e+10, 0.00000e+00, -9.85384e+09, 4.92692e+09,
		-1.95243e+09, -1.26862e+10, 0.00000e+00, -1.02373e+10, 5.11866e+09,
		6.09466e+09, -1.50422e+09, -1.49655e+09, -2.06776e+09, -1.14203e+10,
		0.00000e+00, -1.23173e+10, 0.00000e+00, -2.56008e+09, -1.58681e+09,
		-8.45026e+09, 0.00000e+00, -8.84260e+09, 0.00000e+00, -1.80870e+09,
		-4.11028e+08, -6.57673e+08, -9.37080e+08, -1.03037e+09, -1.10633e+09,
		-2.01594e+09, -2.06769e+09, -2.11599e+09, -2.16157e+09, -1.96997e+09,
		-1.99343e+09, -2.01959e+09, -2.46363e+09, -2.11605e+09, -2.16157e+09,
		-1.61091e+09, -1.63308e+09, -9.48232e+10, -1.34999e+09, -3.53904e+08,
		-7.82931e+09, 0.00000e+00, -5.05490e+09, 5.05490e+09, -1.42645e+09,
		-9.44531e+09, 0.00000e+00, -5.57536e+09, 5.57536e+09, 3.10083e+08,
		1.57470e+09, -1.09239e+09, -1.42212e+09, -7.82931e+09, 0.00000e+00,
		-8.42484e+09, 0.00000e+00, -1.74915e+09, -1.09152e+09, -5.80176e+09,
		0.00000e+00, -6.06332e+09, 0.00000e+00, -1.23945e+09, -2.10657e+08,
		-4.59425e+08, -5.41617e+08, -6.08099e+08, -6.63462e+08, -1.38754e+09,
		-1.42203e+09, -1.45417e+09, -1.48445e+09, -1.41719e+09, -1.41961e+09,
		-1.42795e+09, -1.43955e+09, -1.45424e+09, -1.48451e+09, -1.10759e+09,
		-1.12237e+09, -7.30116e+10, -1.03833e+09, -1.26090e+08, -6.00104e+09,
		-1.28165e+05, -1.21134e+09, -2.57854e+09, 3.86780e+09, -1.18446e+09,
		-7.42533e+09, -7.71512e+08, -3.00807e+09, 4.51210e+09, 1.30296e+08,
		1.88709e+09, 3.37907e+09, -1.09264e+09, -6.00111e+09, 0.00000e+00,
		-6.44634e+09, 0.00000e+00, -1.33727e+09, -8.38734e+08, -4.45178e+09,
		2.38419e-07, -4.64795e+09, 0.00000e+00, -9.49679e+08, -2.21187e+08,
		-2.96109e+08, -3.56846e+08, -4.07369e+08, -1.06663e+09, -1.09253e+09,
		-1.11662e+09, -1.13929e+09, -1.15141e+09, -1.13866e+09, -1.13543e+09,
		-1.13742e+09, -1.11671e+09, -1.13936e+09, -8.50787e+08, -8.61872e+08,
		-8.53257e+09, 8.82980e+08, 4.40658e+09, -3.57443e+09, 3.83154e+09,
		1.61991e+09, 4.47900e+09, 3.39002e+08, 0.00000e+00, 9.69882e+07,
		4.04625e+09, 9.38122e+07, 2.16157e+09, -2.98798e+09, -1.09264e+09,
		-1.44483e+05, -6.05655e+04, -2.62619e+04, -7.45235e+02, 3.05143e+04,
		1.02827e+05, 1.51503e+09, -6.33451e+08, 0.00000e+00, -6.44767e+09,
		5.88402e+08, 5.92193e+09, 6.24161e+09, -3.77543e+09, -7.21818e+07,
		-3.58783e+07, 0.00000e+00, 3.56187e+07, 6.44714e+04, -7.79768e+04,
		-4.68877e+04, 2.69236e+03, 5.09640e+03, -1.39634e+05, -1.10634e+05,
		6.64457e+08, 5.83119e+08, 3.94627e+08, -2.16945e+08, 1.52603e+05,
		9.86642e+04, -9.92470e+09, -9.41595e+08, 5.78674e+09, -2.39864e+09,
		6.54174e+08, 1.93690e+09, 2.11716e+09, 1.00471e+09, 3.49830e+08,
		-2.48254e+08, 4.12768e+05, 3.26260e+05, -1.67704e+04, 1.00577e+05,
		1.19424e+05, 1.11542e+05, 9.19582e+04, 3.78234e+04, -8.06219e+08,
		-7.08071e+08, -6.33451e+08, -2.29870e+08, 1.21903e+04, -1.25085e+10,
		-2.24268e+09, 1.12134e+09, -1.60272e+09, 8.01362e+08, -4.93931e+09,
		-2.56008e+09, -1.74915e+09, -1.33727e+09, 6.40827e+04, 1.28165e+05,
		1.92248e+05, 8.69834e+05, 6.63484e+05, 1.10248e+06, 3.81567e+05,
		1.30493e+05, 0.00000e+00, -8.12176e+04, -6.54749e+09, 7.38242e+08,
		3.41748e+09, -2.58254e+09, 3.02345e+09, 1.28338e+09, 3.55408e+09,
		2.42906e+08, 0.00000e+00, 4.82093e+07, 3.07791e+09, 4.62756e+07,
		1.63308e+09, -2.25188e+09, -8.38734e+08, -1.20194e+05, -4.78960e+04,
		-1.75657e+04, 5.92709e+03, 3.59538e+04, 1.07434e+05, 1.25705e+09,
		-5.09837e+08, 0.00000e+00, -5.00897e+09, 4.76194e+08, 4.49897e+09,
		4.70167e+09, -2.60551e+09, -3.13868e+07, -1.56934e+07, 0.00000e+00,
		1.56934e+07, 5.93884e+04, -6.76680e+04, -3.90856e+04, 7.33255e+03,
		1.01039e+04, -1.22419e+05, -9.67684e+04, 5.36243e+08, 4.68036e+08,
		3.14320e+08, -1.76559e+08, 1.38336e+05, 9.01815e+04, -7.63134e+09,
		-7.71491e+08, 4.21856e+09, -1.89500e+09, 5.16819e+08, 1.37332e+09,
		1.49426e+09, 7.15443e+08, 2.47890e+08, -1.81900e+08, 3.67039e+05,
		2.89487e+05, -8.75075e+03, 9.13216e+04, 1.06865e+05, 9.92939e+04,
		8.15322e+04, 3.39271e+04, -6.56377e+08, -5.73059e+08, -5.09837e+08,
		-1.84060e+08, 1.10737e+04, -8.96744e+09, -1.75740e+09, 8.78701e+08,
		-1.23487e+09, 6.17437e+08, -3.47348e+09, -1.80870e+09, -1.23945e+09,
		-9.49679e+08, 5.71334e+04, 1.14267e+05, 1.71400e+05, 7.49158e+05,
		5.67893e+05, 9.46015e+05, 3.28713e+05, 1.12698e+05, 0.00000e+00,
		-7.03743e+04, -1.08676e+10, -2.33357e+05, 6.82797e+08, -1.69771e+06,
		0.00000e+00, -1.83746e+06, 0.00000e+00, -3.82527e+05, 9.68089e+07,
		-2.66192e+09, 1.95165e+09, -5.46982e+08, 1.09029e+09, -4.32259e+08,
		8.61365e+08, -4.95685e+08, -1.41880e+09, -4.00681e+08, -3.08719e+08,
		-1.87359e+08, 1.43806e+08, 2.69309e+09, 1.22277e+08, 1.45667e+09,
		1.08228e+08, 1.03124e+09, 9.80725e+07, 8.10391e+08, -3.03720e+05,
		5.40668e+08, -1.08473e+09, -1.34006e+09, -2.80155e+05, 4.39896e+08,
		-8.82797e+08, -2.97817e+05, -3.05852e+05, -3.13364e+05, -3.20452e+05,
		-5.01398e+08, -4.39870e+08, -3.93108e+08, -3.11188e+05, -3.19124e+05,
		-2.84695e+05, -2.89335e+05, -1.14201e+10, -1.51778e+05, -8.84453e+05,
		0.00000e+00, -9.53918e+05, 0.00000e+00, -1.98266e+05, -6.17080e+08,
		-2.73862e+09, 9.75827e+08, -5.45908e+08, 5.45526e+08, -4.31340e+08,
		4.31011e+08, -5.46049e+08, -7.09400e+08, 0.00000e+00, 0.00000e+00,
		-1.94826e+08, 1.00916e+08, 2.76937e+09, 9.06005e+07, 1.49655e+09,
		8.29874e+07, 1.05392e+09, 7.70413e+07, 8.25330e+08, -1.60139e+05,
		1.13854e+09, -1.73648e+09, -1.49826e+09, -3.73434e+08, -1.47469e+05,
		9.29408e+08, -1.41899e+09, 5.66405e+05, -1.56132e+05, -1.60139e+05,
		-1.63879e+05, -1.67404e+05, -5.27139e+08, -4.61942e+08, -1.63879e+05,
		-1.67404e+05, -1.49709e+05, -1.51770e+05, -1.17913e+10, -1.04562e+05,
		-6.06344e+05, 0.00000e+00, -6.52465e+05, 0.00000e+00, -1.35464e+05,
		-6.38702e+08, -1.81697e+09, -1.95165e+09, -3.91479e+05, -1.08990e+09,
		-3.38095e+05, -8.61027e+08, -5.81792e+08, 0.00000e+00, 4.00681e+08,
		3.08719e+08, -2.00145e+08, 8.03497e+07, 2.84941e+09, 7.40857e+07,
		1.53052e+09, 6.91443e+07, 1.07374e+09, 1.15769e+09, 1.09253e+09,
		8.38633e+08, 1.80523e+09, -1.92989e+09, -1.69296e+09, -4.25993e+08,
		5.98255e+08, 1.47874e+09, -1.58697e+09, 2.88175e+04, 1.46344e+05,
		-1.07467e+05, -1.10136e+05, -1.12625e+05, -1.14968e+05, -5.45212e+08,
		-1.12625e+05, -1.14968e+05, -1.02934e+05, -1.04307e+05, -1.20975e+10,
		-8.04366e+04, -4.64758e+05, -2.91038e-11, -4.99240e+05, 0.00000e+00,
		-1.03565e+05, -6.55896e+08, 9.14396e+07, -5.85496e+09, 1.09009e+09,
		-3.27057e+09, 8.61192e+08, -2.58384e+09, -6.10388e+08, 7.09400e+08,
		8.01362e+08, 6.17437e+08, -2.04407e+08, 6.78568e+07, 2.91726e+09,
		6.35758e+07, 1.56020e+09, 1.48212e+09, 1.42203e+09, 1.09144e+09,
		1.17375e+09, 1.11662e+09, 8.50708e+08, 2.56289e+09, -1.63908e+09,
		-1.94563e+09, -4.89812e+08, 6.49038e+08, 1.57239e+09, 2.10739e+09,
		-1.36474e+09, 1.21090e+04, 1.75376e+05, 3.14032e+05, -8.26179e+04,
		-8.46197e+04, -8.64840e+04, -8.82386e+04, -8.82386e+04, -8.00977e+04,
		-1.23569e+10, -6.57093e+04, -3.78574e+05, 0.00000e+00, -4.06088e+05,
		0.00000e+00, -8.41840e+04, -2.29162e+09, 2.90359e+09, -1.92642e+09,
		2.40431e+09, -1.61621e+09, 2.98235e+09, -9.75827e+09, 2.72564e+09,
		-5.45129e+09, 2.15334e+09, -4.30668e+09, -1.34401e+09, 0.00000e+00,
		-2.08020e+08, 5.93170e+07, 2.97659e+09, 2.12389e+09, 2.06769e+09,
		1.58674e+09, 1.50771e+09, 1.45417e+09, 1.10753e+09, 1.19054e+09,
		1.13929e+09, 8.61807e+08, -5.71726e+08, 7.15262e+08, 1.72575e+09,
		5.67181e+03, 1.95410e+05, 2.64523e+04, 6.88643e+04, -6.74569e+04,
		-6.90589e+04, -7.05496e+04, -7.19517e+04, -4.03260e+08, -1.25840e+10,
		-5.57477e+04, -3.20414e+05, -1.45519e-11, -3.43293e+05, 0.00000e+00,
		-7.11253e+04, -6.83058e+08, 3.50681e+09, -2.61701e+09, 2.80592e+09,
		-1.96987e+09, 4.92620e+09, -3.95169e+09, 3.84555e+09, -2.23964e+09,
		-5.60669e+08, 3.04424e+09, -1.78166e+09, -4.39350e+08, -1.36522e+09,
		0.00000e+00, -2.11186e+08, 4.00227e+09, 3.94921e+09, 3.02965e+09,
		2.16663e+09, 2.11599e+09, 1.61086e+09, 1.53299e+09, 1.48445e+09,
		1.12232e+09, -2.91610e+09, 8.03718e+08, 1.92165e+09, -1.44519e+09,
		2.60509e+03, 2.10516e+05, 1.53018e+04, 8.07825e+04, 2.23772e+05,
		-5.72062e+04, -5.85420e+04, -5.97841e+04, -6.09516e+04, -4.53940e+08,
		-1.37325e+10, 3.73301e+09, -1.12494e+08, 2.32749e+05, -2.61618e+05,
		-1.75570e+05, -6.74569e+04, 6.87653e+03, -7.94246e+08, -2.30104e+09,
		0.00000e+00, -8.93509e+08, 5.57253e+07, 6.21653e+08, 8.03556e+03,
		3.05852e+05, 7.74679e+03, 1.63879e+05, 7.50103e+03, 1.14968e+05,
		-6.83498e+08, -6.02874e+08, -5.41353e+08, -4.92737e+08, -1.26620e+10,
		1.95980e+09, -1.07957e+08, 3.54018e+05, -2.71407e+05, -2.31405e+05,
		-8.46197e+04, 6.65583e+03, 6.96235e+04, -8.11585e+08, -1.89930e+09,
		9.49650e+08, -1.03758e+08, 3.84550e+09, -5.19604e+08, -1.90125e+08,
		7.51129e+03, 3.13364e+05, 7.26533e+03, 1.67404e+05, -6.98001e+08,
		-6.15354e+08, -5.52318e+08, -5.02524e+08, -1.19740e+10, 1.35014e+09,
		-1.04102e+08, 5.30909e+05, -2.71170e+05, -3.40769e+05, -1.12625e+05,
		6.46766e+03, 8.08799e+04, 1.32464e+05, -8.27656e+08, -1.46466e+09,
		1.46466e+09, -9.69882e+07, 3.94926e+09, -1.00029e+08, 2.01602e+09,
		4.70639e+08, 7.08795e+03, 3.20452e+05, -7.11424e+08, -6.26891e+08,
		-5.62443e+08, -5.11553e+08, -1.17245e+10, 1.03862e+09, -1.00742e+08,
		8.06367e+05, -6.11592e+05, -6.59993e+05, -1.67404e+05, 6.30064e+03,
		9.76647e+04, 1.55192e+05, 6.86587e+04, -1.37317e+05, -1.29610e+05,
		-8.42722e+08, -9.91734e+08, 1.48760e+09, -9.15219e+07, 4.04625e+09,
		-9.38122e+07, 2.06776e+09, -9.68557e+07, 1.38765e+09, -7.23993e+08,
		-6.37682e+08, -5.71905e+08, -5.19983e+08, -4.68370e+10, -6.67579e+08,
		-3.78758e+09, 0.00000e+00, 2.83697e+09, 9.63549e+06, -4.55050e+09,
		-1.12344e+09, -9.31636e+08, -5.71726e+08, -1.96566e+08, -1.40705e+09,
		-5.56030e+08, 4.69040e+07, -3.03006e+09, 1.51503e+09, -4.03110e+09,
		0.00000e+00, -2.51410e+09, 1.25705e+09, -3.28188e+09, 0.00000e+00,
		-6.15285e+07, -7.89173e+08, -1.14731e+08, -7.65315e+08, -1.59706e+08,
		-7.53506e+08, -6.83498e+08, -6.98001e+08, -7.11424e+08, -7.23993e+08,
		-3.78940e+10, -5.89157e+08, -2.66566e+09, 1.33283e+09, 1.27397e+09,
		9.74610e+07, -2.75959e+09, 1.41699e+09, -8.40852e+08, -4.89812e+08,
		-1.20883e+08, 2.62162e+08, -1.06058e+07, 2.98035e+09, -1.00691e+09,
		-7.13606e+08, -2.66566e+09, 1.33283e+09, -3.54036e+09, 0.00000e+00,
		-2.19635e+09, 1.09817e+09, -2.86529e+09, 0.00000e+00, -4.67617e+07,
		-7.20298e+08, -8.85077e+07, -6.98203e+08, -6.02874e+08, -6.15354e+08,
		-6.26891e+08, -6.37682e+08, -3.22241e+10, -5.26710e+08, -1.79118e+09,
		1.79118e+09, 7.73281e+08, 1.69627e+08, -1.18219e+09, 3.11765e+09,
		-7.64966e+08, -4.25993e+08, -6.47721e+07, 3.12760e+08, 7.03327e+08,
		-2.24086e+07, 3.10660e+09, 5.60325e+07, 1.35184e+09, 1.51050e+07,
		-2.38824e+09, 1.19412e+09, -3.16725e+09, 0.00000e+00, -1.95518e+09,
		9.77589e+08, -2.54919e+09, 0.00000e+00, -3.68048e+07, -6.66122e+08,
		-5.41353e+08, -5.52318e+08, -5.62443e+08, -5.71905e+08, -2.58855e+10,
		-4.75189e+08, -1.08473e+09, 1.62709e+09, 5.35387e+08, 2.30028e+08,
		-1.38141e+08, 3.01783e+09, -6.93522e+08, -3.73434e+08, -2.21125e+07,
		3.48956e+08, 7.34812e+08, 2.73044e+08, -5.46088e+08, -6.83410e+08,
		-2.96523e+07, 3.21999e+09, 3.00402e+07, 1.42204e+09, 1.12486e+08,
		8.23387e+08, -2.16945e+09, 1.08473e+09, -2.29870e+09, 1.14935e+09,
		-1.76559e+09, 8.82797e+08, -1.84060e+09, 9.20299e+08, -4.92737e+08,
		-5.02524e+08, -5.11553e+08, -5.19983e+08, -8.73406e+09, 1.76596e+09,
		5.93313e+09, -5.25102e+09, 3.28639e+09, 3.81601e+09, 4.64192e+09,
		1.76642e+09, 5.60669e+08, 9.15219e+07, 4.13777e+09, -4.40012e+09,
		-1.45424e+09, 8.35125e+07, 3.05852e+05, -7.51129e+03, -7.08071e+08,
		0.00000e+00, 5.74676e+08, -6.72091e+04, 3.68854e+04, -3.42574e+04,
		7.81760e+03, 6.10354e+04, -8.92404e+09, 5.60669e+08, 5.60669e+08,
		1.12134e+09, -8.52203e+09, 9.43441e+09, -9.85237e+09, -2.16157e+09,
		8.13559e+07, 1.26108e+09, 1.60139e+05, -7.26533e+03, 2.72242e+09,
		-1.74583e+09, -8.06219e+08, 0.00000e+00, 6.33451e+08, 1.14935e+09,
		-7.08795e+03, 3.13364e+05, -4.01207e+04, 1.38960e+04, -2.04357e+04,
		-4.14882e+03, -6.64794e+09, 1.47648e+09, 4.47900e+09, -3.80222e+09,
		2.59276e+09, 3.01478e+09, 3.69646e+09, 1.35882e+09, 4.39350e+08,
		4.43418e+07, 3.12225e+09, -3.30096e+09, -1.10759e+09, 3.95702e+07,
		2.81564e+05, -4.48031e+03, -5.73059e+08, 0.00000e+00, 4.60149e+08,
		-5.30399e+04, 4.19876e+04, -2.45664e+04, 1.38518e+04, 6.84519e+04,
		-6.74033e+09, 4.39350e+08, 4.39350e+08, 8.78701e+08, -6.35844e+09,
		6.97695e+09, -7.27207e+09, -1.63308e+09, 3.79796e+07, 9.18842e+08,
		1.47469e+05, -4.30060e+03, 2.19942e+09, -1.42566e+09, -6.56377e+08,
		0.00000e+00, 5.09837e+08, 9.20299e+08, -4.12088e+03, 2.86044e+05,
		-2.94596e+04, 1.94790e+04, -1.22948e+04, 2.46324e+03, -1.14200e+10,
		-7.40979e+07, -8.66144e+08, -6.12189e+08, -8.32080e+07, 4.15298e+07,
		-1.34349e+07, 2.64359e+09, -1.03219e+08, 2.06438e+08, -2.58126e+09,
		-2.46662e+09, -2.15230e+09, -1.64895e+09, -3.45131e+05, -2.66725e+05,
		-1.80059e+05, -8.88139e+04, 5.50288e+03, -3.61259e+04, 1.76119e+05,
		-8.77393e+04, 1.42251e+04, -1.24224e+08, -2.00816e+08, -3.43584e+08,
		-7.61231e+08, -7.62588e+09, 4.44286e+08, 4.80128e+08, 1.87251e+08,
		-1.63286e+09, -1.68025e+09, 6.94661e+08, -4.52863e+09, -2.09154e+09,
		-1.13184e+09, -2.75663e+08, 3.09711e+08, -5.58634e+07, -2.26414e+08,
		-2.05274e+08, 4.02451e+09, -2.17488e+05, -2.47270e+05, -2.21838e+05,
		-1.70113e+05, -9.80544e+04, -6.86587e+04, 6.91388e+04, 1.77577e+08,
		-6.57510e+07, -2.30445e+08, -2.55742e+08, -2.09064e+10, -4.44286e+08,
		3.27285e+09, -1.03219e+08, 2.06438e+08, -1.59506e+09, -1.25300e+09,
		-3.79272e+08, -5.21879e+08, -8.93914e+08, -2.63327e+09, -2.77086e+09,
		-2.64724e+09, -2.21791e+09, -2.39864e+09, 6.54174e+08, -2.24268e+08,
		-1.89500e+09, 5.16819e+08, -1.75740e+08, -3.17978e+08, -5.18691e+08,
		-5.66122e+08, -4.36894e+08, -1.10556e+08, -2.88910e+08, -4.60207e+08,
		-1.89930e+08, -9.76440e+07, -9.91734e+07, 1.99185e+08, 2.60935e+08,
		4.22371e+07, -1.03010e+10, -6.24033e+07, 4.25466e+09, 1.00744e+09,
		-1.48495e+09, -1.56187e+09, -4.51673e+09, -2.28794e+09, -1.31316e+09,
		-2.89170e+08, 4.54717e+08, -1.16838e+08, -5.37213e+08, -3.32167e+08,
		-1.60272e+08, -1.23487e+08, -1.55342e+08, -3.32043e+08, -4.39005e+08,
		-4.42773e+08, -3.10405e+08, -2.73044e+08, -5.09020e+06, -1.10842e+10,
		5.88402e+08, -1.58243e+09, 1.88772e+09, -1.59635e+09, -4.45644e+09,
		-3.44404e+09, -2.32497e+09, -1.14679e+09, -2.46350e+04, -1.66686e+05,
		-1.23966e+05, -8.46231e+04, -6.47989e+04, -1.99907e+05, -1.91029e+05,
		-1.41763e+08, -5.49795e+08, -6.82239e+08, -2.16945e+08, -5.76165e+04,
		-2.50264e+04, -1.52271e+10, -9.41595e+08, -2.00908e+09, -2.39864e+09,
		6.54174e+08, -2.80827e+09, 1.11014e+09, -3.19284e+09, -2.86444e+09,
		-2.19656e+09, -1.63130e+05, -8.01711e+04, -3.50721e+05, -1.61980e+05,
		-8.76555e+04, -2.13488e+04, 5.94712e+04, 2.62593e+04, -8.06219e+08,
		-7.08071e+08, -6.33451e+08, -2.29870e+08, -1.21903e+04, -8.38014e+09,
		4.76194e+08, -1.11107e+09, 1.30275e+09, -1.04129e+09, -3.26758e+09,
		-2.54293e+09, -1.73007e+09, -8.66817e+08, -2.26893e+04, -1.49202e+05,
		-1.11865e+05, -7.82916e+04, -5.93307e+04, -1.78334e+05, -1.70770e+05,
		-1.20134e+08, -4.48858e+08, -5.52404e+08, -1.76559e+08, -5.21322e+04,
		-2.28473e+04, -1.14044e+10, -7.71491e+08, -1.40819e+09, -1.89500e+09,
		5.16819e+08, -2.01980e+09, 7.80597e+08, -2.29506e+09, -2.06153e+09,
		-1.58231e+09, -1.44743e+05, -7.05021e+04, -3.06587e+05, -1.42816e+05,
		-7.69442e+04, -1.80074e+04, 5.33825e+04, 2.32063e+04, -6.56377e+08,
		-5.73059e+08, -5.09837e+08, -1.84060e+08, -1.10737e+04};
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
				BOOST_REQUIRE(iter != row.end());
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 20);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 9);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 18);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 1);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 4);
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
	BOOST_REQUIRE_EQUAL(cluster.getId(), 15);
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

BOOST_AUTO_TEST_CASE(traits)
{
	using NetworkType = AlloyReactionNetwork;
	using Spec = NetworkType::Species;

	static_assert(ReactionNetworkTraits<NetworkType>::numSpecies == 6);
	BOOST_REQUIRE_EQUAL(numberOfSpecies<Spec>(), 6);
	BOOST_REQUIRE_EQUAL(numberOfInterstitialSpecies<Spec>(), 3);
	BOOST_REQUIRE_EQUAL(numberOfVacancySpecies<Spec>(), 3);

	// toLabelString
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::V), "V");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::Void), "Void");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::Faulted), "Faulted");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::I), "I");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::Perfect), "Perfect");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::Frank), "Frank");

	// toNameString
	BOOST_REQUIRE_EQUAL(toNameString(Spec::V), "Vacancy");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::Void), "Void");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::Faulted), "Faulted");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::I), "Interstitial");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::Perfect), "Perfect");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::Frank), "Frank");

	// parseSpeciesId / getSpeciesLabel / getSpeciesName
    auto network = NetworkType();
	auto sid = network.parseSpeciesId("V");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::V);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "V");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Vacancy");
	sid = network.parseSpeciesId("Void");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::Void);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "Void");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Void");
	sid = network.parseSpeciesId("Faulted");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::Faulted);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "Faulted");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Faulted");
	sid = network.parseSpeciesId("I");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::I);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "I");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Interstitial");
	sid = network.parseSpeciesId("Perfect");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::Perfect);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "Perfect");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Perfect");
	sid = network.parseSpeciesId("Frank");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::Frank);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "Frank");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Frank");

	// Define the sequence
	using GroupingRange = SpeciesForGrouping<Spec, 6>;

	// It should always return 0 because they are all orthogonal
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::V), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::Void), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::Faulted), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::I), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::Perfect), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::Frank), 0);
}

BOOST_AUTO_TEST_SUITE_END()
