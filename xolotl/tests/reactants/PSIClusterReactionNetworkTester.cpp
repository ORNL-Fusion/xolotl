#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <PSIHeVCluster.h>
#include <PSIHeCluster.h>
#include <PSIVCluster.h>
#include <PSIInterstitialCluster.h>
#include <PSIHeInterstitialCluster.h>
#include <xolotlPerf.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the ReactionNetwork
 */
BOOST_AUTO_TEST_SUITE(PSIReactionNetwork_testSuite)

BOOST_AUTO_TEST_CASE(checkReactants) {
	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>(registry);

	// Add a few He, V and I PSIClusters
	auto heCluster = std::unique_ptr<PSIHeCluster>(
			new PSIHeCluster(10, *(psiNetwork.get()), registry));
	auto vCluster = std::unique_ptr<PSIVCluster>(
			new PSIVCluster(4, *(psiNetwork.get()), registry));
	auto interstitialCluster = std::unique_ptr<PSIInterstitialCluster>(
			new PSIInterstitialCluster(48, *(psiNetwork.get()), registry));
	psiNetwork->add(std::move(heCluster));
	psiNetwork->add(std::move(vCluster));
	psiNetwork->add(std::move(interstitialCluster));

	// Check the network, He first
	auto retHeCluster = (PSICluster *) psiNetwork->get(Species::He, 10);
	BOOST_REQUIRE(retHeCluster);
	BOOST_REQUIRE_EQUAL("He_10", retHeCluster->getName());
	BOOST_REQUIRE_EQUAL(10, retHeCluster->getSize());
	// V
	auto retVCluster = (PSICluster *) psiNetwork->get(Species::V, 4);
	BOOST_REQUIRE(retVCluster);
	BOOST_REQUIRE_EQUAL(4, retVCluster->getSize());
	BOOST_REQUIRE_EQUAL("V_4", retVCluster->getName());
	// I
	auto retICluster = (PSICluster *) psiNetwork->get(Species::I, 48);
	BOOST_REQUIRE(retICluster);
	BOOST_REQUIRE_EQUAL(48, retICluster->getSize());
	BOOST_REQUIRE_EQUAL("I_48", retICluster->getName());

	// Check the getter for all reactants
	auto& clusters = psiNetwork->getAll();
	BOOST_REQUIRE_EQUAL(3U, clusters.size());
	// Check the size of the network
	BOOST_REQUIRE_EQUAL(3, psiNetwork->size());

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
	int maxClusterSize = 10;
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			auto cluster = std::unique_ptr<PSIHeVCluster>(
					new PSIHeVCluster(numHe, numV, *(psiNetwork.get()),
							registry));
			psiNetwork->add(std::move(cluster));
			counter++;
		}
	}

	BOOST_TEST_MESSAGE("Added " << counter << " HeV clusters");

	// Add a whole bunch of HeI clusters to make sure that the network can
	// handle large numbers of them properly too. Use a different max cluster
	// size to throw a little variability at the network.
	counter = 0;
	maxClusterSize = 9;
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			auto cluster = std::unique_ptr<PSIHeInterstitialCluster>(
					new PSIHeInterstitialCluster(numHe, numI,
							*(psiNetwork.get()), registry));
			psiNetwork->add(std::move(cluster));
			counter++;
		}
	}

	BOOST_TEST_MESSAGE("Added " << counter << " HeI clusters");

	// Try adding a duplicate HeV and catch the exception
	auto duplicateCluster = std::unique_ptr<PSIHeVCluster>(
			new PSIHeVCluster(5, 3, *(psiNetwork.get()), registry));
	try {
		psiNetwork->add(std::move(duplicateCluster));
		BOOST_FAIL(
				"Test failed because adding a duplicate" << " to the network was allowed.");
	} catch (const std::string& /* e */) {
		// Do nothing. It was supposed to fail.
	}

	// Make sure that everything was added
	auto& reactants = psiNetwork->getAll();
	BOOST_REQUIRE_EQUAL(84U, reactants.size());

	// Get the clusters by type and check them. Start with He.
	auto& heReactants = psiNetwork->getAll(ReactantType::He);
	auto& heReactant = heReactants.begin()->second;
	BOOST_REQUIRE_EQUAL(1U, heReactants.size());
	BOOST_REQUIRE_EQUAL("He_10", heReactant->getName());
	// V
	auto& vReactants = psiNetwork->getAll(ReactantType::V);
	auto& vReactant = vReactants.begin()->second;
	BOOST_REQUIRE_EQUAL(1U, vReactants.size());
	BOOST_REQUIRE_EQUAL("V_4", vReactant->getName());
	// I
	auto& iReactants = psiNetwork->getAll(ReactantType::I);
	auto& iReactant = iReactants.begin()->second;
	BOOST_REQUIRE_EQUAL(1U, iReactants.size());
	BOOST_REQUIRE_EQUAL("I_48", iReactant->getName());

	// HeV
	auto& heVReactants = psiNetwork->getAll(ReactantType::HeV);
	BOOST_REQUIRE_EQUAL(45U, heVReactants.size());
	// HeI
	auto& heIReactants = psiNetwork->getAll(ReactantType::HeI);
	BOOST_REQUIRE_EQUAL(36U, heIReactants.size());

	// Add the required He_1, V_1, I_1 clusters to the network.
	heCluster = std::unique_ptr<PSIHeCluster>(
			new PSIHeCluster(1, *(psiNetwork.get()), registry));
	vCluster = std::unique_ptr<PSIVCluster>(
			new PSIVCluster(1, *(psiNetwork.get()), registry));
	interstitialCluster = std::unique_ptr<PSIInterstitialCluster>(
			new PSIInterstitialCluster(1, *(psiNetwork.get()), registry));
	psiNetwork->add(std::move(heCluster));
	psiNetwork->add(std::move(vCluster));
	psiNetwork->add(std::move(interstitialCluster));

	// Try changing the temperature and make sure it works
	psiNetwork->setTemperature(1000.0);
	IReactant& reactant = reactants.at(0);
	BOOST_REQUIRE_CLOSE(1000.0, reactant.getTemperature(), 0.0001);

	return;
}

BOOST_AUTO_TEST_CASE(checkProperties) {
	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>(registry);

	// Check the properties
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::He), 0);
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::HeV), 0);
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::HeI), 0);

	// Add a couple of clusters
	auto heCluster = std::unique_ptr<PSIHeCluster>(
			new PSIHeCluster(5, *(psiNetwork.get()), registry));
	psiNetwork->add(std::move(heCluster));
	auto heVCluster = std::unique_ptr<PSIHeVCluster>(
			new PSIHeVCluster(5, 3, *(psiNetwork.get()), registry));
	psiNetwork->add(std::move(heVCluster));

	// Check the properties again
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::He), 5);
	BOOST_REQUIRE_EQUAL(psiNetwork->getMaxClusterSize(ReactantType::HeV), 8);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {
	// Create the network
	auto psiNetwork = make_shared<PSIClusterReactionNetwork>(registry);

	// Check the names of all cluster types. Use a simple counting
	// system to look over the list.
	auto names = psiNetwork->getKnownReactantTypes();
	unsigned int marker = 0;
	for (auto name : names) {
		if (name == ReactantType::He)
			++marker;
		else if (name == ReactantType::D)
			++marker;
		else if (name == ReactantType::T)
			++marker;
		else if (name == ReactantType::V)
			++marker;
		else if (name == ReactantType::I)
			++marker;
		else if (name == ReactantType::PSIMixed)
			++marker;
		else if (name == ReactantType::PSISuper)
			++marker;
		else if (name == ReactantType::HeV)
			++marker;
		else if (name == ReactantType::HeI)
			++marker;
	}
	BOOST_REQUIRE_EQUAL(9U, marker);
	BOOST_REQUIRE_EQUAL(marker, names.size());

	return;
}

/**
 * This operation tests the operations of the ReactionNetwork that copy the
 * concentrations to and from a client array.
 */
BOOST_AUTO_TEST_CASE(checkArrayOperations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork();
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
