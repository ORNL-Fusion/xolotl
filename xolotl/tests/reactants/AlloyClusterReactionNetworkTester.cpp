#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <AlloyCluster.h>
#include "SimpleReactionNetwork.h"
#include <AlloyVacCluster.h>
#include <AlloyIntCluster.h>
#include <xolotlPerf/xolotlPerf.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the ReactionNetwork
 */
BOOST_AUTO_TEST_SUITE(AlloyReactionNetwork_testSuite)

BOOST_AUTO_TEST_CASE(checkReactants) {
	// Create the network
	auto alloyNetwork = make_shared<AlloyClusterReactionNetwork>(registry);

	// Add a vacancy AlloyCluster
	auto vCluster = std::unique_ptr<AlloyVacCluster>(
			new AlloyVacCluster(5, *(alloyNetwork.get()), registry));
	alloyNetwork->add(std::move(vCluster));

	// Check the network
	auto retVCluster = (AlloyCluster *) alloyNetwork->get(Species::V, 5);
	BOOST_REQUIRE(retVCluster);
	BOOST_REQUIRE_EQUAL("V_5", retVCluster->getName());
	BOOST_REQUIRE_EQUAL(5, retVCluster->getSize());

	// Check the getter for all reactants
	auto& clusters = alloyNetwork->getAll();
	BOOST_REQUIRE_EQUAL(1U, clusters.size());
	// Check the size of the network
	BOOST_REQUIRE_EQUAL(1, alloyNetwork->size());

	// Check the cluster ids.
	int id = retVCluster->getId();
	BOOST_REQUIRE(id > 0 && id <= 1);

	// Try adding a duplicate V and catch the exception
	auto duplicateCluster = std::unique_ptr<AlloyVacCluster>(
			new AlloyVacCluster(5, *(alloyNetwork.get()), registry));
	try {
		alloyNetwork->add(std::move(duplicateCluster));
		BOOST_FAIL(
				"Test failed because adding a duplicate to the network was allowed.");
	} catch (const std::string& /* e */) {
		// Do nothing. It was supposed to fail.
	}

	// Make sure that everything was added
	auto& reactants = alloyNetwork->getAll();
	BOOST_REQUIRE_EQUAL(1U, reactants.size());

	// Try changing the temperature and make sure it works
	//Set a diffusion factor so that the temperature matters
	retVCluster->setDiffusionFactor(1.0);
	// Add a grid point for the temperature
	alloyNetwork->addGridPoints(1);
	alloyNetwork->setTemperature(1000.0, 0);
	BOOST_REQUIRE_CLOSE(1000.0, retVCluster->getTemperature(0), 0.0001);

	return;
}

BOOST_AUTO_TEST_CASE(checkProperties) {
	// Create the network
	auto alloyNetwork = make_shared<AlloyClusterReactionNetwork>(registry);

	// Check the properties
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Void), 0);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Frank), 0);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Faulted), 0);
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::Perfect), 0);

	// Add a cluster
	auto iCluster = std::unique_ptr<AlloyIntCluster>(
			new AlloyIntCluster(3, *(alloyNetwork.get()), registry));
	alloyNetwork->add(std::move(iCluster));

	// Check the properties again
	BOOST_REQUIRE_EQUAL(alloyNetwork->getMaxClusterSize(ReactantType::I), 3);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {
	// Create the network
	auto alloyNetwork = make_shared<AlloyClusterReactionNetwork>(registry);

	// Check the names of all cluster types. Use a simple counting
	// system to look over the list.
	auto names = alloyNetwork->getKnownReactantTypes();
	unsigned int marker = 0;
	for (auto name : names) {
		if (name == ReactantType::V)
			++marker;
		else if (name == ReactantType::I)
			++marker;
		else if (name == ReactantType::Void)
			++marker;
		else if (name == ReactantType::Faulted)
			++marker;
		else if (name == ReactantType::Perfect)
			++marker;
		else if (name == ReactantType::Frank)
			++marker;
		else if (name == ReactantType::VoidSuper)
			++marker;
		else if (name == ReactantType::FaultedSuper)
			++marker;
		else if (name == ReactantType::PerfectSuper)
			++marker;
		else if (name == ReactantType::FrankSuper)
			++marker;
	}
	BOOST_REQUIRE_EQUAL(10U, marker);
	BOOST_REQUIRE_EQUAL(marker, names.size());

	return;
}

/**
 * This operation tests the operations of the ReactionNetwork that copy the
 * concentrations to and from a client array.
 */
BOOST_AUTO_TEST_CASE(checkArrayOperations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();
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
