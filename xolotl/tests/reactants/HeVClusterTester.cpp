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
#include <HeCluster.h>
#include <HeVCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <math.h>
#include <limits>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry = std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the HeVCluster.
 */
BOOST_AUTO_TEST_SUITE(HeVCluster_testSuite)


BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	HeVCluster cluster(4, 5, registry);

	// Get the composition back
	auto composition = cluster.getComposition();

	// Check the composition is the created one
	BOOST_REQUIRE_EQUAL(composition["He"], 4);
	BOOST_REQUIRE_EQUAL(composition["V"], 5);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);
}

/**
 * This operation checks the ability of the HeVCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	auto props = network->getProperties();
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the HeV cluster
	// with 3He and 2V
	
	{
		// Get the connectivity array from the reactant
		vector<int> composition = {3, 2, 0 };
		auto reactant = (PSICluster *) network->getCompound("HeV", composition);
		// Check the type name
		BOOST_REQUIRE_EQUAL("HeV",reactant->getType());
		auto reactionConnectivity = reactant->getConnectivity();
		
		BOOST_REQUIRE_EQUAL(reactant->getComposition().at("He"), 3);
		BOOST_REQUIRE_EQUAL(reactant->getComposition().at("V"), 2);
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			
			// V
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			
			// I
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			
			// HeV
			0, 0, 1, 0, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0,
			0, 0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 1,
			0, 0,
			0,
			
			// HeI
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0,
			0, 0,
			0
		};
		
		for (int i = 0; i < reactionConnectivity.size(); i++) {
			BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
		}
	}
}

/**
 * This operation checks the ability of the HeVCluster to compute the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

	// Get an HeV cluster with compostion 2,1,0.
	vector<int> composition = {2, 1, 0};
	auto cluster = (PSICluster *) network->getCompound(
			"HeV",composition);
	// Get one that it combines with (He)
	auto secondCluster = (PSICluster *) network->get("He", 1);
	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(0.0);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
	vector<double> energies = {3.02, 7.25,
			numeric_limits<double>::infinity(), 10.2};
	cluster->setTemperature(1000.0);
	cluster->setBindingEnergies(energies);
	cluster->setConcentration(0.5);

	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(2.950E+10);
	secondCluster->setMigrationEnergy(0.13);
	secondCluster->setTemperature(1000.0);
	energies = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(),
			numeric_limits<double>::infinity(), 8.27};
	secondCluster->setBindingEnergies(energies);
	secondCluster->setConcentration(0.5);

 	// Compute the rate constants that are needed for the flux
 	cluster->computeRateConstants(1000.0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(1000.0);
	BOOST_TEST_MESSAGE("HeVClusterTester Message: \n" << "Total Flux is " << flux << "\n"
			  << "   -Production Flux: " << cluster->getProductionFlux(1000.0) << "\n"
			  << "   -Combination Flux: " << cluster->getCombinationFlux(1000.0) << "\n"
			  << "   -Dissociation Flux: " << cluster->getDissociationFlux(1000.0) << "\n"
	  	  	  << "   -Emission Flux: " << cluster->getEmissionFlux(1000.0) << "\n");
	BOOST_REQUIRE_CLOSE(-8964899015.0, flux, 10.0);
}

/**
 * This operation checks the reaction radius for HeVCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	vector<shared_ptr<HeVCluster>> clusters;
	shared_ptr<HeVCluster> cluster;
	double expectedRadii[] = { 0.1372650265, 0.1778340462, 0.2062922619,
			0.2289478080, 0.2480795532 };

	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<HeVCluster>(new HeVCluster(1, i, registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				.000001);
	}
}

BOOST_AUTO_TEST_SUITE_END()
