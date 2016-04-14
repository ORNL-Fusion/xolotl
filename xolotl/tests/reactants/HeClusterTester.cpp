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
 * This suite is responsible for testing the HeCluster.
 */
BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/**
 * This operation checks the ability of the HeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	auto props = network->getProperties();
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the 6th He reactant (numHe=6)
	// Get the connectivity array from the reactant
	auto reactant = (PSICluster *) network->get("He", 6);

	// Check the type name
	BOOST_REQUIRE_EQUAL("He",reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();
		
	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
		// He
		1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
			
		// V
		1, 1, 1, 1, 0, 0, 0, 0, 0, 0,

		// I
		1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
		// HeV
		1, 1, 1, 0, 0, 0, 0, 0, 0,
		1, 1, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0,
		0, 0,
		0,
			
		// HeI
		1, 1, 1, 0, 0, 0, 0, 0, 0,
		1, 1, 0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0,
		0, 0,
		0
	};

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the HeCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

	// Get an He cluster with composition 1,0,0.
	auto cluster = (PSICluster *) network->get("He", 1);
	// Get one that it combines with (He2)
	auto secondCluster = (PSICluster *) network->get("He", 2);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.950E+10);
	cluster->setMigrationEnergy(0.13);
 	cluster->setTemperature(1000.0);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(3.240E+010);
	secondCluster->setMigrationEnergy(0.2);
	secondCluster->setConcentration(0.5);
 	secondCluster->setTemperature(1000.0);

 	// Compute the rate constants that are needed for the flux
 	cluster->computeRateConstants();
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux();
	BOOST_TEST_MESSAGE("HeClusterTester Message: \n" << "Total Flux is " << flux << "\n"
			  << "   -Production Flux: " << cluster->getProductionFlux() << "\n"
			  << "   -Combination Flux: " << cluster->getCombinationFlux() << "\n"
			  << "   -Dissociation Flux: " << cluster->getDissociationFlux() << "\n"
			  << "   -Emission Flux: " << cluster->getEmissionFlux() << "\n");

	BOOST_REQUIRE_CLOSE(6110430723517.8, flux, 0.1);

	return;
}

/**
 * This operation checks the HeCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = {-1.96821e+11, 1.23317e+13, 3.21149e+12,
			-1.79298e+10, -1.95933e+10, 0.0, -1.87741e+10, -2.04376e+10, 0.0, 2.23350e+12,
			2.25143e+12, 2.46031e+12, -1.79298e+10, 2.25143e+12, 0.0};
	// Get the simple reaction network
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork(3);

	// Get an He cluster with composition 1,0,0.
	auto cluster = (PSICluster *) network->get("He", 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.950E+10);
	cluster->setMigrationEnergy(0.13);
 	cluster->setTemperature(1000.0);
	cluster->setConcentration(0.5);

 	// Compute the rate constants that are needed for the partial derivatives
 	cluster->computeRateConstants();
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives();

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 15U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	return;
}

/**
 * This operation checks the reaction radius for HeCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create a helium cluster
	shared_ptr<HeCluster> cluster;

	// The vector of radii to compare with
	double expectedRadii[] = { 0.3, 0.3237249066, 0.3403673722, 0.3536164159,
				0.3648047284, 0.3745846085, 0.3833299460, 0.3912773576,
				0.3985871973, 0.4053737480 };

	// Check all the values
	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<HeCluster>(new HeCluster(i, registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i-1], cluster->getReactionRadius(), 0.000001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
