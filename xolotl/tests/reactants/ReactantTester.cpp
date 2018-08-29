#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <Reactant.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the Reactant.
 */BOOST_AUTO_TEST_SUITE(Reactant_testSuite)

BOOST_AUTO_TEST_CASE(checkTemperature) {
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);
	// Create a reactant
	Reactant reactant(*(network.get()), registry);
	// Add rates
	reactant.addGridPoints(1);
	// Set its temperature
	reactant.setTemperature(1000.0, 0);

	// Check its temperature
	BOOST_REQUIRE_CLOSE(1000.0, reactant.getTemperature(0), 0.0001);

	return;
}

BOOST_AUTO_TEST_CASE(checkComposition) {
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);
	// Add rates
	network->addGridPoints(1);
	// Create a reactant
	Reactant reactant(*(network.get()), registry);

	// Check its default composition
	BOOST_REQUIRE_EQUAL(10U, reactant.getComposition().size());

	return;
}

BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Create a network
	auto network = make_shared<PSIClusterReactionNetwork>(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Add rates
	network->addGridPoints(1);
	// Create a reactant
	Reactant reactant(*(network.get()), registry);

	// Check its default partial derivatives
	BOOST_REQUIRE_EQUAL(1U, reactant.getPartialDerivatives(0).size());

	// Create a reference and temp partial derivative vector
	std::vector<double> refPartials = std::vector<double>(3, 0.0);
	std::vector<double> tempPartials = std::vector<double>(3, 0.0);

	reactant.getPartialDerivatives(tempPartials, 0);

	// Compare to the refence one
	BOOST_REQUIRE_EQUAL(tempPartials[0], refPartials[0]);
	BOOST_REQUIRE_EQUAL(tempPartials[1], refPartials[1]);
	BOOST_REQUIRE_EQUAL(tempPartials[2], refPartials[2]);

	return;
}

BOOST_AUTO_TEST_CASE(checkCopying) {
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);
	// Create a reference reactant
	Reactant reactant(*(network.get()), registry);
	// Add rates
	reactant.addGridPoints(1);
	reactant.setId(5);
	reactant.setConcentration(10.0);
	reactant.setTemperature(5.0, 0);

	// Copy the Reactant
	Reactant reactantCopy(reactant);

	// Check the ids and names
	BOOST_REQUIRE_EQUAL(reactant.getId(), reactantCopy.getId());
	BOOST_REQUIRE_EQUAL(reactant.getName(), reactantCopy.getName());
	BOOST_REQUIRE_CLOSE(5.0, reactant.getTemperature(0), 0.0001);

	// Increase the concentration
	reactantCopy.setConcentration(15.0);

	// The values should now be different,
	// so check them against the known values
	BOOST_REQUIRE_CLOSE(reactant.getConcentration(), 10.0, 1.0e-7);
	BOOST_REQUIRE_CLOSE(reactantCopy.getConcentration(), 15.0, 1.0e-7);

	return;
}

BOOST_AUTO_TEST_CASE(checkConcentration) {
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);
	// Add rates
	network->addGridPoints(1);
	// Create a reactant
	Reactant reactant(*(network.get()), registry);
	reactant.setConcentration(1.0);

	// Make sure it was set correctly
	BOOST_REQUIRE_EQUAL(1.0, reactant.getConcentration());

	// Make sure the base class getTotalFlux returns 0 for now
	BOOST_REQUIRE_EQUAL(0.0, reactant.getTotalFlux(0));

	return;
}

BOOST_AUTO_TEST_SUITE_END()
