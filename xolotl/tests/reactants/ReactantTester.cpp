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
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the Reactant.
 */
BOOST_AUTO_TEST_SUITE(Reactant_testSuite)

/**
 * This operation tests the copy constructor.
 */
BOOST_AUTO_TEST_CASE(checkCopying) {
	
	// Create a reference Reactant
	std::shared_ptr<Reactant> reactant(new Reactant);
	reactant->setConcentration(10.0);
	
	// Copy the Reactant
	std::shared_ptr<Reactant> reactant2(new Reactant(*reactant));
	
	// Check that the pointers are different
	BOOST_REQUIRE_NE(reactant.get(), reactant2.get());
	
	reactant2->increaseConcentration(5.0);
	
	// The values should now be different,
	// so check them against the known values
	BOOST_REQUIRE_CLOSE(reactant->getConcentration(), 10.0, 1e-7);
	BOOST_REQUIRE_CLOSE(reactant2->getConcentration(), 15.0, 1e-7);
}

BOOST_AUTO_TEST_CASE(checkManipulateConcentration) {

	// Create a Reactant
	std::shared_ptr<Reactant> reactant(new Reactant);
	reactant->setConcentration(1.0);

	// Make sure it was set correctly
	BOOST_CHECK_EQUAL(1.0, reactant->getConcentration());

	// Increase it
	reactant->increaseConcentration(3.3);

	// Make sure its correct
	BOOST_CHECK_EQUAL(4.3, reactant->getConcentration());

	// Decrease it
	reactant->decreaseConcentration(1.3);

	// Make sure its correct
	BOOST_CHECK_EQUAL(3.0, reactant->getConcentration());

	// Zero it
	reactant->zero();

	// Check it was zeroed
	BOOST_CHECK_EQUAL(0.0, reactant->getConcentration());

	// Make sure the base class getTotalFlux returns 0 for now
	BOOST_CHECK_EQUAL(0.0, reactant->getTotalFlux(0.0));

}


/*BOOST_AUTO_TEST_CASE(toClusterMap) {
	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	
	std::map<std::string, int> cluster;
	
	// Test a couple of the HeClusters
	
	cluster = network->toClusterMap(0);
	BOOST_REQUIRE_EQUAL(cluster["He"], 1);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	cluster = network->toClusterMap(9);
	BOOST_REQUIRE_EQUAL(cluster["He"], 10);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	// Test VClusters
	
	cluster = network->toClusterMap(10);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 1);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	cluster = network->toClusterMap(19);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 10);
	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
	
	// Test IClusters
	
	cluster = network->toClusterMap(20);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 1);
	
	cluster = network->toClusterMap(29);
	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
	BOOST_REQUIRE_EQUAL(cluster["I"], 10);
	
	// Test HeVClusters
	
	for (int i = 30; i < 75; i++) {
		// Get the actual He and V amounts
		
		shared_ptr<Reactant> reactant = network->reactants->at(i);
		shared_ptr<HeVCluster> cluster = std::dynamic_pointer_cast<HeVCluster>(reactant);
		int actualHe = cluster->getSpeciesSize("He");
		int actualV = cluster->getSpeciesSize("V");
		
		// Get the amounts expected by Reactant::toClusterMap()
		
		std::map<std::string, int> species = network->toClusterMap(i);
		int expectedHe = species["He"];
		int expectedV = species["V"];
		
		BOOST_REQUIRE_EQUAL(actualHe, expectedHe);
		BOOST_REQUIRE_EQUAL(actualV, expectedV);
	}
}*/

BOOST_AUTO_TEST_CASE(checkIsConnected) {

	//shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	// shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();

	/*for (int i = 0; i < network->reactants->size(); i++) {
		for (int j = 0; j < network->reactants->size(); j++) {
			//std::cout << i << " " << j << "\n";
			//BOOST_REQUIRE_EQUAL(network->isConnected(i,j), network->isConnected(j,i));
		}
	}*/
}
BOOST_AUTO_TEST_SUITE_END()

