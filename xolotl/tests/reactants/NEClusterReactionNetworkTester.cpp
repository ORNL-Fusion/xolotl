#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <NECluster.h>
#include "SimpleReactionNetwork.h"
#include <NEXeCluster.h>
#include <xolotlPerf.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the ReactionNetwork
 */
BOOST_AUTO_TEST_SUITE(NEReactionNetwork_testSuite)

BOOST_AUTO_TEST_CASE(checkReactants) {
	// Create the network
	auto neNetwork = make_shared<NEClusterReactionNetwork>(registry);

	// Add a Xe NECluster
	auto xeCluster = std::unique_ptr<NEXeCluster>(
			new NEXeCluster(10, *(neNetwork.get()), registry));
	neNetwork->add(std::move(xeCluster));

	// Check the network, Xe first
	auto retXeCluster = (NECluster *) neNetwork->get(Species::Xe, 10);
	BOOST_REQUIRE(retXeCluster);
	BOOST_REQUIRE_EQUAL("Xe_10", retXeCluster->getName());
	BOOST_REQUIRE_EQUAL(10, retXeCluster->getSize());

	// Check the getter for all reactants
	auto& clusters = neNetwork->getAll();
	BOOST_REQUIRE_EQUAL(1U, clusters.size());
	// Check the size of the network
	BOOST_REQUIRE_EQUAL(1, neNetwork->size());

	// Check the cluster ids.
	int id = retXeCluster->getId();
	BOOST_REQUIRE(id > 0 && id <= 1);

	// Try adding a duplicate Xe and catch the exception
	auto duplicateCluster = std::unique_ptr<NEXeCluster>(
			new NEXeCluster(10, *(neNetwork.get()), registry));
	try {
		neNetwork->add(std::move(duplicateCluster));
		BOOST_FAIL(
				"Test failed because adding a duplicate" << " to the network was allowed.");
	} catch (const std::string& /* e */) {
		// Do nothing. It was supposed to fail.
	}

	// Make sure that everything was added
	auto& reactants = neNetwork->getAll();
	BOOST_REQUIRE_EQUAL(1U, reactants.size());

	// Try changing the temperature and make sure it works
	neNetwork->setTemperature(1000.0);
	IReactant& reactant = reactants.at(0);
	BOOST_REQUIRE_CLOSE(1000.0, reactant.getTemperature(), 0.0001);

	return;
}

BOOST_AUTO_TEST_CASE(checkProperties) {
	// Create the network
	auto neNetwork = make_shared<NEClusterReactionNetwork>(registry);

	// Check the properties
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::Xe), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::V), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::I), 0);
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::XeV), 0);

	// Add a cluster
	auto xeCluster = std::unique_ptr<NEXeCluster>(
			new NEXeCluster(5, *(neNetwork.get()), registry));
	neNetwork->add(std::move(xeCluster));

	// Check the properties again
	BOOST_REQUIRE_EQUAL(neNetwork->getMaxClusterSize(ReactantType::Xe), 5);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {
	// Create the network
	auto neNetwork = make_shared<NEClusterReactionNetwork>(registry);

	// Check the names of all cluster types. Use a simple counting
	// system to look over the list.
	auto names = neNetwork->getKnownReactantTypes();
	unsigned int marker = 0;
	for (auto name : names) {
		if (name == ReactantType::Xe)
			++marker;
		else if (name == ReactantType::V)
			++marker;
		else if (name == ReactantType::I)
			++marker;
		else if (name == ReactantType::NESuper)
			++marker;
		else if (name == ReactantType::XeV)
			++marker;
		else if (name == ReactantType::XeI)
			++marker;
	}
	BOOST_REQUIRE_EQUAL(6U, marker);
	BOOST_REQUIRE_EQUAL(marker, names.size());

	return;
}

/**
 * This operation tests the operations of the ReactionNetwork that copy the
 * concentrations to and from a client array.
 */
BOOST_AUTO_TEST_CASE(checkArrayOperations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleNEReactionNetwork();
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
