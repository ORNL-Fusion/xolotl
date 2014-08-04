/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <HeInterstitialCluster.h>
#include <HeCluster.h>
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
 * This suite is responsible for testing the HeInterstitialCluster.
 */
BOOST_AUTO_TEST_SUITE(HeInterstitialCluster_testSuite)

BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	HeInterstitialCluster cluster(4, 2, registry);

	// Get the composition back
	auto composition = cluster.getComposition();

	// Check the composition is the created one
	BOOST_REQUIRE_EQUAL(composition["He"], 4);
	BOOST_REQUIRE_EQUAL(composition["V"], 0);
	BOOST_REQUIRE_EQUAL(composition["I"], 2);
}

/**
 * This operation checks the ability of the HeInterstitialCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	auto props = network->getProperties();

	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";

	// Check the reaction connectivity of the HeI cluster
	// with 5He and 3I

	{
		// Get the connectivity array from the reactant
		vector<int> composition = { 5, 0, 3 };
		auto reactant =
				(PSICluster *) (network->getCompound("HeI", composition));
		// Check the type name
		BOOST_REQUIRE_EQUAL("HeI", reactant->getType());
		auto reactionConnectivity = reactant->getConnectivity();

		BOOST_REQUIRE_EQUAL(reactant->getComposition().at("He"), 5);
		BOOST_REQUIRE_EQUAL(reactant->getComposition().at("I"), 3);

		// Check the connectivity for He, V, and I

		int connectivityExpected[] = {
				// He
				1, 1, 1, 1, 1, 0, 0, 0, 0, 0,

				// V
				1, 1, 0, 0, 0, 0, 0, 0, 0, 0,

				// I
				1, 0, 1, 0, 0, 0, 0, 0, 0, 0,

				// HeV
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0,

				// HeI
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1,
				1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0 };

		for (int i = 0; i < reactionConnectivity.size(); i++) {
			BOOST_REQUIRE_EQUAL(reactionConnectivity[i],
					connectivityExpected[i]);
		}
	}
}

/**
 * This operation checks the ability of the HeInterstitialCluster to compute
 * the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {
	BOOST_TEST_MESSAGE(
			"HeInterstitialClusterTester Message: \n" << "BOOST_AUTO_TEST_CASE(checkTotalFlux): " << "Arbitrary values used because of lack of data!" << "\n");

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

	// Get an HeI cluster with compostion 1,0,1.
	vector<int> composition = { 1, 0, 1 };
	auto cluster = (PSICluster *) network->getCompound("HeI", composition);
	// Get one that it combines with (I)
	auto secondCluster = (PSICluster *) network->get("I", 1);
	// Set the diffusion factor, migration and binding energies to arbitrary
	// values because HeI does not exist in benchmarks
	cluster->setDiffusionFactor(1.5E+10);
	cluster->setTemperature(1000.0);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
	vector<double> energies = { 5.09, numeric_limits<double>::infinity(), 5.09,
			12.6 };
	cluster->setBindingEnergies(energies);
	cluster->setConcentration(0.5);

	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(2.13E+10);
	secondCluster->setMigrationEnergy(0.013);
	energies = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(),
		numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
	secondCluster->setBindingEnergies(energies);
	secondCluster->setConcentration(0.5);
	secondCluster->setTemperature(1000.0);

	// Compute the rate constants that are needed for the flux
	cluster->computeRateConstants(1000.0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(1000.0);
	BOOST_TEST_MESSAGE(
			"HeInterstitialClusterTester Message: \n" << "Total Flux is " << flux << "\n" << "   -Production Flux: " << cluster->getProductionFlux(1000.0) << "\n" << "   -Combination Flux: " << cluster->getCombinationFlux(1000.0) << "\n" << "   -Dissociation Flux: " << cluster->getDissociationFlux(1000.0) << "\n" << "   -Emission Flux: " << cluster->getEmissionFlux(1000.0) << "\n");

	// Check the flux
	BOOST_REQUIRE_CLOSE(-16982855380.0, flux, 10.0);

	return;
}

/**
 * This operation checks the HeCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0};
	// Get the simple reaction network
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork(3);

	// Get an HeI cluster with compostion 2,0,1.
	vector<int> composition = { 2, 0, 1 };
	auto cluster = (PSICluster *) network->getCompound("HeI", composition);
	// Set the diffusion factor, migration and binding energies to arbitrary
	// values because HeI does not exist in benchmarks
	cluster->setDiffusionFactor(1.5E+10);
	cluster->setTemperature(1000.0);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
	vector<double> energies = { 5.09, numeric_limits<double>::infinity(), 5.09,
			12.6 };
	cluster->setBindingEnergies(energies);
	cluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the partial derivatives
	cluster->computeRateConstants(1000.0);
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(1000.0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 15);

	// Check all the values
	for (int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 10.0);
	}
}

/**
 * This operation checks the reaction radius for HeInterstitialCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	vector<shared_ptr<HeInterstitialCluster>> clusters;
	shared_ptr<HeInterstitialCluster> cluster;
	double expectedRadii[] = { 0.1372650265, 0.1778340462, 0.2062922619,
			0.2289478080, 0.2480795532 };

	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr < HeInterstitialCluster
				> (new HeInterstitialCluster(1, i, registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				.000001);
	}
}

BOOST_AUTO_TEST_SUITE_END()
