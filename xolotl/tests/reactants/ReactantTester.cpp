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

BOOST_AUTO_TEST_SUITE_END()

