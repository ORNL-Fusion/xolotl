#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <NECluster.h>
#include "SimpleReactionNetwork.h"
#include <XeCluster.h>
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

	// Add a few Xe NEClusters
	auto xeCluster = make_shared<XeCluster>(10, registry);
	neNetwork->add(xeCluster);

	// Check the network, Xe first
	auto retXeCluster = (NECluster *) neNetwork->get("Xe", 10);
	BOOST_REQUIRE(retXeCluster);
	BOOST_REQUIRE_EQUAL("Xe_10", retXeCluster->getName());
	BOOST_REQUIRE_EQUAL(10, retXeCluster->getSize());

	// Check the getter for all reactants
	auto clusters = neNetwork->getAll();
	BOOST_REQUIRE_EQUAL(1U, clusters->size());
	// Check the size of the network
	BOOST_REQUIRE_EQUAL(1, neNetwork->size());

	// Check the cluster ids.
	int id = retXeCluster->getId();
	BOOST_REQUIRE(id > 0 && id <= 1);

	// Try adding a duplicate Xe and catch the exception
	shared_ptr<XeCluster> duplicateCluster = std::make_shared<XeCluster>(10,
			registry);
	try {
		neNetwork->add(duplicateCluster);
		BOOST_FAIL(
				"Test failed because adding a duplicate" << " to the network was allowed.");
	} catch (const std::string& /* e */) {
		// Do nothing. It was supposed to fail.
	}

	// Make sure that everything was added
	auto reactants = neNetwork->getAll();
	BOOST_REQUIRE_EQUAL(1U, reactants->size());

	// Try changing the temperature and make sure it works
	neNetwork->setTemperature(1000.0);
	BOOST_REQUIRE_CLOSE(1000.0, reactants->at(0)->getTemperature(), 0.0001);

	return;
}

BOOST_AUTO_TEST_CASE(checkProperties) {
	// Create the network
	auto neNetwork = make_shared<NEClusterReactionNetwork>(registry);

	// Access network "properties."
	auto numXeClusters = neNetwork->getNumXeClusters();
	auto numVClusters = neNetwork->getNumVClusters();
	auto numIClusters = neNetwork->getNumIClusters();
	auto numXeVClusters = neNetwork->getNumXeVClusters();
	auto numXeIClusters = neNetwork->getNumXeIClusters();
	auto maxXeVClusterSize = neNetwork->getMaxXeVClusterSize();
	auto maxXeIClusterSize = neNetwork->getMaxXeIClusterSize();
	auto maxXeClusterSize = neNetwork->getMaxXeClusterSize();
	auto maxVClusterSize = neNetwork->getMaxVClusterSize();
	auto maxIClusterSize = neNetwork->getMaxIClusterSize();

	// Check the properties
	BOOST_REQUIRE_EQUAL(0, numXeClusters);
	BOOST_REQUIRE_EQUAL(0, numVClusters);
	BOOST_REQUIRE_EQUAL(0, numIClusters);
	BOOST_REQUIRE_EQUAL(0, numXeVClusters);
	BOOST_REQUIRE_EQUAL(0, numXeIClusters);
	BOOST_REQUIRE_EQUAL(0, maxXeVClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxXeIClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxXeClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxVClusterSize);
	BOOST_REQUIRE_EQUAL(0, maxIClusterSize);

	// Add a cluster
	auto heCluster = make_shared<XeCluster>(5, registry);
	neNetwork->add(heCluster);

	// Grab the properties afresh
	numXeClusters = neNetwork->getNumXeClusters();
	maxXeClusterSize = neNetwork->getMaxXeClusterSize();

	// Check the properties again
	BOOST_REQUIRE_EQUAL(1, numXeClusters);
	BOOST_REQUIRE_EQUAL(5, maxXeClusterSize);

	return;
}

BOOST_AUTO_TEST_CASE(checkNames) {
	// Create the network
	auto neNetwork = make_shared<NEClusterReactionNetwork>(registry);

	// Check the names of the regular cluster types. Use a simple counting
	// system to look over the list since there is no way to check exact
	// containment with a vector.
	auto names = neNetwork->getNames();
	unsigned int marker = 0;
	for (unsigned int i = 0; i < names.size(); i++) {
		if (names[i] == "Xe")
			++marker;
		else if (names[i] == "V")
			++marker;
		else if (names[i] == "I")
			++marker;
		else if (names[i] == "NESuper")
			++marker;
	}
	BOOST_REQUIRE_EQUAL(4U, marker);
	BOOST_REQUIRE_EQUAL(marker, names.size());

	// Check the names of the compound cluster types. Use the same counting
	// system as above.
	auto compoundNames = neNetwork->getCompoundNames();
	marker = 0;
	for (unsigned int i = 0; i < compoundNames.size(); i++) {
		if (compoundNames[i] == "XeV")
			++marker;
		else if (compoundNames[i] == "XeI")
			++marker;
	}
	BOOST_REQUIRE_EQUAL(2U, marker);
	BOOST_REQUIRE_EQUAL(marker, compoundNames.size());

	return;
}

/**
 * This operation tests the copy constructor.
 */
BOOST_AUTO_TEST_CASE(checkCopying) {
	//NEClusterReactionNetwork network;
	NEClusterReactionNetwork network(registry);

	// Add a reactant
	shared_ptr<Reactant> heCluster(new XeCluster(1, registry));
	heCluster->setConcentration(50.0);
	network.add(heCluster);

	// Copy the network
	NEClusterReactionNetwork networkCopy = network;

	// Check that the ReactionNetwork fields are copied
	BOOST_REQUIRE_EQUAL(network.getNumXeClusters(),
			networkCopy.getNumXeClusters());

	// Check that changing the concentration of a copy does not update the
	// original. Start by updating the copy.
	auto copiedXeCluster = networkCopy.get("Xe", 1);
	copiedXeCluster->setConcentration(7.0);
	BOOST_REQUIRE_CLOSE(7.0, copiedXeCluster->getConcentration(), 1.0e-5);
	// Make sure the original wasn't changed.
	BOOST_REQUIRE_CLOSE(50.0, heCluster->getConcentration(), 1.0e-5);

	// Check the size of the network
	BOOST_REQUIRE_EQUAL(1, networkCopy.size());

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
		BOOST_REQUIRE_CLOSE(1.0, reactants->at(0)->getConcentration(), 1.0e-15);
	}

	// Clear memory
	delete[] concentrations;

	return;
}

BOOST_AUTO_TEST_CASE(checkRefCounts) {
	// Obtain a network to work with.
	// This network was built programmatically.
	shared_ptr<ReactionNetwork> network = getSimpleNEReactionNetwork();

	// Because each Reactant in the network is given a pointer
	// (a shared_ptr) to the network, and we have one shared_ptr to it,
	// its reference count should be network's size + 1.
	BOOST_TEST_MESSAGE("After creation, network size: " << network->size());
	BOOST_TEST_MESSAGE(
			"After creation, network ref count: " << network.use_count());
	BOOST_REQUIRE_EQUAL(network.use_count(), network->size() + 1);

	// Tell the network to break dependency cycles between
	// the Reactants in the network and the network itself.
	// In a "real" use, this allows the network and Reactants
	// to be destroyed gracefully when the shared_ptr pointing
	// to the network goes out of scope, because it allows
	// the network's reference count to reach zero.
	network->askReactantsToReleaseNetwork();

	// All objects from within the network should have released their
	// shared_ptr to the network, so our shared_ptr should be the
	// only remaining shared_ptr.  Thus, the network's reference
	// count should be 1 at this point.
	// If it is, when our shared_ptr goes out of scope the network will
	// be destroyed.  We can't easily test that it is destroyed.
	BOOST_TEST_MESSAGE(
			"After releasing network refs, network size: " << network->size());
	BOOST_TEST_MESSAGE(
			"After releasing network refs, network ref count: " << network.use_count());
	BOOST_REQUIRE_EQUAL(network.use_count(), 1);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
