#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <FeCluster.h>
#include "SimpleReactionNetwork.h"
#include <FeHeVCluster.h>
#include <FeHeCluster.h>
#include <FeVCluster.h>
#include <FeInterstitialCluster.h>
#include <xolotlPerf/xolotlPerf.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the ReactionNetwork
 */
BOOST_AUTO_TEST_SUITE(FeReactionNetwork_testSuite)

BOOST_AUTO_TEST_CASE(checkReactants) {
	// Create the network
	auto feNetwork = make_shared<FeClusterReactionNetwork>(registry);

	// Add a few He, V and I FeClusters
	auto heCluster = std::unique_ptr<FeHeCluster>(
			new FeHeCluster(10, *(feNetwork.get()), registry));
	auto vCluster = std::unique_ptr<FeVCluster>(
			new FeVCluster(4, *(feNetwork.get()), registry));
	auto interstitialCluster = std::unique_ptr<FeInterstitialCluster>(
			new FeInterstitialCluster(1, *(feNetwork.get()), registry));
	feNetwork->add(std::move(heCluster));
	feNetwork->add(std::move(vCluster));
	feNetwork->add(std::move(interstitialCluster));

	// Check the network, He first
	auto retHeCluster = (FeCluster *) feNetwork->get(Species::He, 10);
	BOOST_REQUIRE(retHeCluster);
	BOOST_REQUIRE_EQUAL("He_10", retHeCluster->getName());
	BOOST_REQUIRE_EQUAL(10, retHeCluster->getSize());
	// V
	auto retVCluster = (FeCluster *) feNetwork->get(Species::V, 4);
	BOOST_REQUIRE(retVCluster);
	BOOST_REQUIRE_EQUAL(4, retVCluster->getSize());
	BOOST_REQUIRE_EQUAL("V_4", retVCluster->getName());
	// I
	auto retICluster = (FeCluster *) feNetwork->get(Species::I, 1);
	BOOST_REQUIRE(retICluster);
	BOOST_REQUIRE_EQUAL(1, retICluster->getSize());
	BOOST_REQUIRE_EQUAL("I_1", retICluster->getName());

	// Check the getter for all reactants
	auto& clusters = feNetwork->getAll();
	BOOST_REQUIRE_EQUAL(3U, clusters.size());
	// Check the size of the network
	BOOST_REQUIRE_EQUAL(3, feNetwork->size());

	// Check the cluster ids. All we can do is check that they are between 1
	// and 5. Start with the He cluster
	int id = retHeCluster->getId();
	BOOST_REQUIRE(id > 0 && id <= 5);
	// V
	id = retVCluster->getId();
	BOOST_REQUIRE(id > 0 && id <= 5);
	// I
	id = retICluster->getId();
	BOOST_REQUIRE(id > 0 && id <= 5);

	// Add a whole bunch of HeV clusters to make sure that the network can
	// handle large numbers of them properly.
	int counter = 0;
	int maxClusterSize = 6;
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
			auto cluster = std::unique_ptr<FeHeVCluster>(
					new FeHeVCluster(numHe, numV, *(feNetwork.get()),
							registry));
			feNetwork->add(std::move(cluster));
			counter++;
		}
	}

	// Try adding a duplicate HeV and catch the exception
	auto duplicateCluster = std::unique_ptr<FeHeVCluster>(
			new FeHeVCluster(5, 3, *(feNetwork.get()), registry));
	try {
		feNetwork->add(std::move(duplicateCluster));
		BOOST_FAIL(
				"Test failed because adding a duplicate" << " to the network was allowed.");
	} catch (const std::string& /* e */) {
		// Do nothing. It was supposed to fail.
	}

	// Make sure that everything was added
	auto& reactants = feNetwork->getAll();
	BOOST_REQUIRE_EQUAL(39U, reactants.size());

	// Get the clusters by type and check them. Start with He.
	auto& heReactants = feNetwork->getAll(ReactantType::He);
	auto& heReactant = heReactants.begin()->second;
	BOOST_REQUIRE_EQUAL(1U, heReactants.size());
	BOOST_REQUIRE_EQUAL("He_10", heReactant->getName());
	// V
	auto& vReactants = feNetwork->getAll(ReactantType::V);
	auto& vReactant = vReactants.begin()->second;
	BOOST_REQUIRE_EQUAL(1U, vReactants.size());
	BOOST_REQUIRE_EQUAL("V_4", vReactant->getName());
	// I
	auto& iReactants = feNetwork->getAll(ReactantType::I);
	auto& iReactant = iReactants.begin()->second;
	BOOST_REQUIRE_EQUAL(1U, iReactants.size());
	BOOST_REQUIRE_EQUAL("I_1", iReactant->getName());

	// HeV
	auto& heVReactants = feNetwork->getAll(ReactantType::HeV);
	BOOST_REQUIRE_EQUAL(36U, heVReactants.size());

	// Add the required He_1, V_1 clusters to the network.
	heCluster = std::unique_ptr<FeHeCluster>(
			new FeHeCluster(1, *(feNetwork.get()), registry));
	vCluster = std::unique_ptr<FeVCluster>(
			new FeVCluster(1, *(feNetwork.get()), registry));
	feNetwork->add(std::move(heCluster));
	feNetwork->add(std::move(vCluster));

	// Try changing the temperature and make sure it works
	// Set the diffusion factor because temperature is only relevant for mobile clusters
	heReactant->setDiffusionFactor(1.0);
	// Add a grid point for the temperature
	feNetwork->addGridPoints(1);
	feNetwork->setTemperature(1000.0, 0);
	BOOST_REQUIRE_CLOSE(1000.0, heReactant->getTemperature(0), 0.0001);

	return;
}

BOOST_AUTO_TEST_CASE(checkProperties) {
	// Create the network
	auto feNetwork = make_shared<FeClusterReactionNetwork>(registry);

	// Check the properties
	BOOST_REQUIRE_EQUAL(feNetwork->getMaxClusterSize(ReactantType::He), 0);
	BOOST_REQUIRE_EQUAL(feNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(feNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(feNetwork->getMaxClusterSize(ReactantType::HeV), 0);

	// Add a couple of clusters
	auto heCluster = std::unique_ptr<FeHeCluster>(
			new FeHeCluster(5, *(feNetwork.get()), registry));
	feNetwork->add(std::move(heCluster));
	auto heVCluster = std::unique_ptr<FeHeVCluster>(
			new FeHeVCluster(5, 3, *(feNetwork.get()), registry));
	feNetwork->add(std::move(heVCluster));

	// Check the properties again
	BOOST_REQUIRE_EQUAL(feNetwork->getMaxClusterSize(ReactantType::He), 5);
	BOOST_REQUIRE_EQUAL(feNetwork->getMaxClusterSize(ReactantType::HeV), 8);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {
	// Create the network
	auto feNetwork = make_shared<FeClusterReactionNetwork>(registry);

	// Check the names of all cluster types. Use a simple counting
	// system to look over the list.
	auto names = feNetwork->getKnownReactantTypes();
	unsigned int marker = 0;
	for (auto name : names) {
		if (name == ReactantType::He)
			++marker;
		else if (name == ReactantType::V)
			++marker;
		else if (name == ReactantType::I)
			++marker;
		else if (name == ReactantType::FeSuper)
			++marker;
		else if (name == ReactantType::HeV)
			++marker;
	}
	BOOST_REQUIRE_EQUAL(5U, marker);
	BOOST_REQUIRE_EQUAL(marker, names.size());

	return;
}

/**
 * This operation tests the operations of the ReactionNetwork that copy the
 * concentrations to and from a client array.
 */
BOOST_AUTO_TEST_CASE(checkArrayOperations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleFeReactionNetwork();
	int size = network->size();
	double * concentrations = new double[size];

	// Set default values
	for (int i = 0; i < size; i++) {
		concentrations[i] = 1.0;
	}

	// Set the array to the values in the Reaction network, which should
	// all be zero.
	network->fillConcentrationsArray(concentrations);
	for (int i = 0; i < size; i++) {
		BOOST_REQUIRE_CLOSE(0.0, concentrations[i], 1.0e-15);
	}

	// Reset values to something else
	for (int i = 0; i < size; i++) {
		concentrations[i] = 1.0;
	}

	// Update the network and check it
	network->updateConcentrationsFromArray(concentrations);
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		IReactant& reactant = reactants.at(i);
		BOOST_REQUIRE_CLOSE(1.0, reactant.getConcentration(), 1.0e-15);
	}

	// Clear memory
	delete[] concentrations;

	return;
}

BOOST_AUTO_TEST_SUITE_END()
