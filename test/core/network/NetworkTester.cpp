#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/ConfOptions.h>

using namespace std;
using namespace xolotl::core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the networks.
 */
BOOST_AUTO_TEST_SUITE(Network_testSuite)

using network_types =
#if defined(__clang__)
	boost::mpl::list<
#elif defined(__GNUC__) || defined(__GNUG__)
	std::tuple<
#endif
		NEReactionNetwork, AlloyReactionNetwork, FeReactionNetwork,
		PSIReactionNetwork<PSIFullSpeciesList>>;

BOOST_AUTO_TEST_CASE_TEMPLATE(defaut_constructor, T, network_types)
{
	// Create a file with reaction data for NE.
	std::ofstream reactionFile("reaction.dat");
	reactionFile
		<< "0 0 1 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "0 1 0 -0.160789024638881 -12.7387253625273 4.2488E+11 4.23317"
		<< std::endl
		<< "0 2 0 2.87060634573234 -8.37649260139129 9.9337E+11 3.43830674"
		<< std::endl
		<< "1 0 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 1 0 3.14296517286615 -1.4324963007117 0 0" << std::endl
		<< "1 2 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 3 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 4 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 5 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 6 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 7 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "1 8 0 11.2375337876032 22.3998792567472 4.0767E+11 4.08453089"
		<< std::endl
		<< "Reactions" << std::endl
		<< "0 0 1 0.00596" << std::endl
		<< "0 1 0 0.00596" << std::endl
		<< "0 2 0 0.00596" << std::endl;
	reactionFile.close();

	using NetworkType = T;
	constexpr auto speciesRange = NetworkType::getSpeciesRange();
	std::vector<typename NetworkType::AmountType> maxSpeciesAmounts;

	for (auto i : speciesRange) {
		maxSpeciesAmounts.push_back(1);
	}
	NetworkType network =
		NetworkType(maxSpeciesAmounts, 1, xolotl::options::ConfOptions{});

	// Remove the file because it is not needed anymore
	std::remove("reaction.dat");

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
