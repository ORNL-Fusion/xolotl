/*
 * PSIClusterNetworkLoaderTester.cpp
 *
 *  Created on: Mar 30, 2013
 *      Author: jaybilly
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the PSIClusterNetworkLoader. It
 * creates a string stream that contains each of the available PSICluster types
 * and checks that the loader returns a list with each type in it.
 */
BOOST_AUTO_TEST_SUITE (PSIClusterNetworkLoader_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkLoading) {

	// Local Declarations
	shared_ptr<stringstream> networkStream(
			new stringstream(stringstream::in | stringstream::out));
	string singleHeString = "1 0 0 0.0 Infinity Infinity 8.2699999999999996\n";
	string singleVString =
			"0 50 0 Infinity 2.4900000000000002 Infinity Infinity\n";
	string singleIString = "0 0 1 Infinity Infinity Infinity Infinity\n";
	string mixedString =
			"1 50 0 6.1600000000000001 2.4900000000000002 Infinity Infinity\n";
	PSIClusterNetworkLoader loader = PSIClusterNetworkLoader();
	shared_ptr<ReactionNetwork> network;

	// Load the network stream. This simulates a file with single He, single
	// V, single I and one mixed-species cluster.
	*networkStream << singleHeString << singleVString << singleIString << mixedString;

	// Diagnostic information
	BOOST_TEST_MESSAGE("CLUSTER DATA");
	BOOST_TEST_MESSAGE("He: " << singleHeString);
	BOOST_TEST_MESSAGE("V: " << singleVString);
	BOOST_TEST_MESSAGE("I: " << singleIString);
	BOOST_TEST_MESSAGE("Mixed: " << mixedString);
	BOOST_TEST_MESSAGE("Full Network data: \n" << (*networkStream).str());

	// Setup the Loader
	loader.setInputstream(networkStream);

	// Load the network

	// FIXME - Does this work? shared_ptr != ptr !!!

	network = loader.load();

	// Check the network. It should not be empty
	BOOST_REQUIRE(!(network->properties->empty()));
	BOOST_REQUIRE(!(network->reactants->empty()));
	// It should have four reactants
	BOOST_REQUIRE_EQUAL(4,network->reactants->size());
	// It should have six properties
	BOOST_REQUIRE_EQUAL(6,network->properties->size());

	// Check the reactants

	// Check the properties

	BOOST_FAIL("Not yet implemented.");

}

BOOST_AUTO_TEST_SUITE_END()
