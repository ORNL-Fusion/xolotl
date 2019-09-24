#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <NECluster.h>
#include "SimpleReactionNetwork.h"
#include <NEXeCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the NEXeCluster.
 */
BOOST_AUTO_TEST_SUITE (NEXeCluster_testSuite)

/**
 * This operation checks the ability of the NEXeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimpleNEReactionNetwork();

	// Check the reaction connectivity of the 6th Xe reactant (numXe=6)
	// Get the connectivity array from the reactant
	auto reactant = (NECluster *) network->get(Species::Xe, 6);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Xe == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Xe
	int connectivityExpected[] = {
	// Xe
			1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the NEXeCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	auto network = getSimpleNEReactionNetwork();
	// Set a fission rate for the diffusion to work
	network->setFissionRate(8.0e-9);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an Xe cluster with compostion 1,0,0.
	auto cluster = (NECluster *) network->get(Species::Xe, 1);
	// Get one that it combines with (Xe2)
	auto secondCluster = (NECluster *) network->get(Species::Xe, 2);
	// Set the diffusion factor and migration energy
	cluster->setDiffusionFactor(2.950E+10);
	cluster->setMigrationEnergy(0.13);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy
	secondCluster->setDiffusionFactor(3.240E+010);
	secondCluster->setMigrationEnergy(0.2);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_TEST_MESSAGE(
			"NEXeClusterTester Message: \n" << "Total Flux is " << flux << "\n" << "   -Production Flux: " << cluster->getProductionFlux(0) << "\n" << "   -Combination Flux: " << cluster->getCombinationFlux(0) << "\n" << "   -Dissociation Flux: " << cluster->getDissociationFlux(0) << "\n" << "   -Emission Flux: " << cluster->getEmissionFlux(0) << "\n");

	BOOST_REQUIRE_CLOSE(1.1555563138, flux, 0.0001);

	return;
}

/**
 * This operation checks the NEXeCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = {-0.235828, 2.46149, 0.686708, 0};
	// Get the simple reaction network
	auto network = getSimpleNEReactionNetwork(3);
	// Set a fission rate for the diffusion to work
	network->setFissionRate(8.0e-9);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an Xe cluster with compostion 1,0,0.
	auto cluster = (NECluster *) network->get(Species::Xe, 1);
	// Set the diffusion factor and migration energy
	cluster->setDiffusionFactor(2.950E+10);
	cluster->setMigrationEnergy(0.13);
	cluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the partial derivatives
	network->setTemperature(1000.0, 0);
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 4U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.001);
	}

	return;
}

/**
 * This operation checks the reaction radius for NEXeCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create a helium cluster
	shared_ptr<NEXeCluster> cluster;

	// Get the simple reaction network
	auto network = getSimpleNEReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.3, 0.36083562, 0.413053674, 0.4546243938,
			0.48972928, 0.5204150186, 0.5478547329, 0.57279084, 0.5957264838,
			0.6170202317 };

	// Check all the values
	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<NEXeCluster>(
				new NEXeCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.000001);
	}

	return;
}

/**
 * This operation checks the reaction radius for NEXeCluster with a different density.
 */
BOOST_AUTO_TEST_CASE(checkDensity) {
	// Create a helium cluster
	shared_ptr<NEXeCluster> cluster;

	// Get the simple reaction network with a density of 10.0 Xe nm-3
	std::shared_ptr<xolotlPerf::IHandlerRegistry> registry = std::make_shared<
			xolotlPerf::DummyHandlerRegistry>();
	auto network = getSimpleNEReactionNetwork(0, registry, 10.0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.3, 0.362783, 0.415283, 0.457078, 0.492373,
			0.523224, 0.550812, 0.575882, 0.598942, 0.62035 };

	// Check all the values
	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<NEXeCluster>(
				new NEXeCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.0001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
