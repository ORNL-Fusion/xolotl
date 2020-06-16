#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/PSIReactionNetwork.h>

using namespace std;
using namespace xolotl::core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the networks.
 */
BOOST_AUTO_TEST_SUITE(Network_testSuite)

using network_types = std::tuple<NEReactionNetwork, AlloyReactionNetwork,
		FeReactionNetwork, PSIReactionNetwork<PSIFullSpeciesList>>;

BOOST_AUTO_TEST_CASE_TEMPLATE( defaut_constructor, T, network_types )
{
	using NetworkType = T;
	constexpr auto speciesRange = NetworkType::getSpeciesRange();
	std::vector<typename NetworkType::AmountType> maxSpeciesAmounts;

	for (auto i : speciesRange) {
		maxSpeciesAmounts.push_back(1);
	}
	NetworkType network = NetworkType(maxSpeciesAmounts, 1,
			xolotl::options::Options { });

	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);

	BOOST_REQUIRE(network.getNumClusters() > 0);
	BOOST_REQUIRE(network.getDOF() > 0);
	BOOST_REQUIRE(network.getDeviceMemorySize() > 0);

	double latticeParameter = 0.5;
	network.setLatticeParameter(latticeParameter);
	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), latticeParameter, 0.01);

	BOOST_REQUIRE(network.getAtomicVolume() > 0.0);

	double impurityRadius = 0.3;
	network.setImpurityRadius(impurityRadius);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), impurityRadius, 0.01);

	double interstitialBias = 1.2;
	network.setInterstitialBias(interstitialBias);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), interstitialBias, 0.01);

	double fissionRate = 1.0e-9;
	network.setFissionRate(fissionRate);
	BOOST_REQUIRE_CLOSE(network.getFissionRate(), fissionRate, 0.01);

	bool enableStdReaction = true;
	network.setEnableStdReaction(enableStdReaction);
	BOOST_REQUIRE(network.getEnableStdReaction() == enableStdReaction);

	bool enableReSolution = true;
	network.setEnableReSolution(enableReSolution);
	BOOST_REQUIRE(network.getEnableReSolution() == enableReSolution);

	typename NetworkType::IndexType gridSize = 10;
	network.setGridSize(gridSize);
	BOOST_REQUIRE(network.getGridSize() == gridSize);

	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE(bounds.size() > 0);

	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE(phaseSpace.size() > 0);

	BOOST_REQUIRE(network.getNumberOfSpecies() > 0);
	BOOST_REQUIRE(network.getNumberOfSpeciesNoI() > 0);

	// TODO: get the subpaving? Get the concentrations?
}

BOOST_AUTO_TEST_SUITE_END()
