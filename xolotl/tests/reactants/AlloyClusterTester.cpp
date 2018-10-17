#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>
#include <AlloyCluster.h>
#include <AlloyClusterReactionNetwork.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>
#include "SimpleReactionNetwork.h"

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the AlloyCluster.
 */
BOOST_AUTO_TEST_SUITE (AlloyCluster_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkDiffusionCoefficient) {
	// Get the simple reaction network
	auto network = getSimpleAlloyReactionNetwork(0);
	// Create a cluster
	AlloyCluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature
	cluster.addGridPoints(1);

	// Check E_m = 0.0
	cluster.setMigrationEnergy(0.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), exp(0.0), 0.00001);
	BOOST_REQUIRE_CLOSE(1.0, cluster.getTemperature(0), 0.0001);

	// Make sure the diffusion coefficient is 0.0 if E_m is infinite
	cluster.setMigrationEnergy(numeric_limits<double>::infinity());
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0, 0.000001);

	// Make sure the diffusion coefficient is zero if the diffusion factor is zero
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(0.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0, 0.000001);

	// Make sure the diffusion coefficient is equal to the diffusion factor
	// if the temperature is infinite
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(numeric_limits<double>::infinity(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 1.0, 0.000001);

	// Throw something random in there to be certain
	cluster.setMigrationEnergy(0.013);
	cluster.setDiffusionFactor(1.08E10);
	cluster.setTemperature(1500.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 9766651101.800613,
			0.0000001);

	return;
}

/**
 * This operation tests the default values returned by select flux routines.
 */
BOOST_AUTO_TEST_CASE(checkDefaultFluxes) {
	// Get the simple reaction network
	auto network = getSimpleAlloyReactionNetwork(0);
	// Create a cluster
	AlloyCluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature
	cluster.addGridPoints(1);

	// Check the default values of the fluxes
	BOOST_REQUIRE_CLOSE(cluster.getProductionFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getCombinationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getDissociationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getTotalFlux(0), 0.0, 1e-5);

	return;
}

/**
 * This operation checks the different properties of the vacancy cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyVacancy) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first V reactant (numV=1)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::V, 1);

	// Check the type name
	BOOST_REQUIRE(ReactantType::V == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for V
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (V2)
	auto secondReactant = (AlloyCluster *) network->get(Species::V, 2);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-1.1037070506226931e-05, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -3.65158e-05, -7.63251e-06, -7.9213e-06,
			-8.1512e-06, -8.34535e-06, 0.000694247, -0.012871, -0.013358,
			-0.0137457, -4.25753e-06, -4.33341e-06, -4.40236e-06, -4.46578e-06,
			-4.52466e-06, -4.57974e-06, -4.63157e-06, -4.6806e-06, -4.72717e-06,
			-4.77158e-06, -4.81405e-06, -4.85479e-06, -4.89396e-06, -4.9317e-06,
			0.0, -4.66135e-06, -4.809e-06, -4.94644e-06, -5.07553e-06,
			-5.19762e-06, -5.31374e-06, -5.42469e-06, -5.53112e-06,
			-5.63352e-06, -5.73232e-06, -5.82788e-06, -5.9205e-06, -6.01043e-06,
			-6.0979e-06, 0, -3.15972e+08, -2.72116e+08, -2.40183e+08,
			-2.15813e+08, -1.96553e+08, -1.80915e+08, -1.67942e+08,
			-1.56989e+08, -1.47606e+08, -1.3947e+08, -1.32339e+08, -1.26033e+08,
			-1.20412e+08, -1.15366e+08, -1.10809e+08, -1.0667e+08, -4.50078e-06,
			-4.66135e-06, -4.809e-06, -4.94644e-06, -5.07553e-06, -5.19762e-06,
			-5.31374e-06, -5.42469e-06, -5.53112e-06, -5.63352e-06,
			-5.73232e-06, -5.82788e-06, -5.9205e-06, -6.01043e-06, -6.0979e-06,
			-6.1831e-06, 0.0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of V1
	BOOST_REQUIRE_CLOSE(0.14069, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the interstitial cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyInterstitial) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first I reactant (numI=1)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::I, 1);

	// Check the type name
	BOOST_REQUIRE(ReactantType::I == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for I
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (I2)
	auto secondReactant = (AlloyCluster *) network->get(Species::I, 2);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-0.037213, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { 0.000694247, -0.012871, -0.013358, -0.0137457,
			-0.0140731, -0.123119, -0.0257343, -0.026708, -0.0274832, -0.014355,
			-0.0146108, -0.0148433, -0.0150571, -0.0152557, -0.0154414,
			-0.0156161, -0.0157814, -0.0159385, -0.0160882, -0.0162314,
			-0.0163687, -0.0165008, -0.0166281, -0.0167509, -0.0157165,
			-0.0162144, -0.0166778, -0.017113, -0.0175246, -0.0179162,
			-0.0182903, -0.0186491, -0.0189944, -0.0193275, -0.0196497,
			-0.019962, -0.0202652, -0.0205601, -0.0208474, -3.15972e+08,
			-2.72116e+08, -2.40183e+08, -2.15813e+08, -1.96553e+08,
			-1.80915e+08, -1.67942e+08, -1.56989e+08, -1.47606e+08, -1.3947e+08,
			-1.32339e+08, -1.26033e+08, -1.20412e+08, -1.15366e+08,
			-1.10809e+08, 0, -0.0151751, -0.0157165, -0.0162144, -0.0166778,
			-0.017113, -0.0175246, -0.0179162, -0.0182903, -0.0186491,
			-0.0189944, -0.0193275, -0.0196497, -0.019962, -0.0202652,
			-0.0205601, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of I1
	BOOST_REQUIRE_CLOSE(0.14069, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the void cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyVoid) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first Void reactant (numV=6)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Void, 6);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Void == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Void
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (V1)
	auto secondReactant = (AlloyCluster *) network->get(Species::V, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-2.128764e-06, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -4.25753e-06, -4.46337e-06, -4.60777e-06,
			-4.72272e-06, 3.52556e-06, -0.014355, -0.015049, -0.0155359,
			-0.0159235, -4.25753e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.63101e+08,
			-3.11389e+08, -2.73846e+08, -2.45268e+08, -2.22735e+08,
			-2.04479e+08, -1.89364e+08, -1.76626e+08, -1.65733e+08,
			-1.56301e+08, -1.48049e+08, -1.40761e+08, -1.34273e+08,
			-1.28458e+08, -1.23211e+08, -1.18453e+08, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Void6
	BOOST_REQUIRE_CLOSE(0.2556446, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the faulted cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyFaulted) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first Faulted reactant (numV=6)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Faulted, 6);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Faulted == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Faulted
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (V1)
	auto secondReactant = (AlloyCluster *) network->get(Species::V, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-2.3306727e-06, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -4.66135e-06, -4.86719e-06, -5.01159e-06,
			-5.12654e-06, -5.22361e-06, -0.0157165, -0.0164106, -0.0168974,
			-0.017285, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			-4.66135e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.9251e+08,
			-3.35897e+08, -2.94853e+08, -2.63649e+08, -2.39074e+08,
			-2.19184e+08, -2.02731e+08, -1.88879e+08, -1.77044e+08,
			-1.66805e+08, -1.57852e+08, -1.49951e+08, -1.42923e+08,
			-1.36627e+08, -1.30951e+08, -1.25805e+08, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Faulted6
	BOOST_REQUIRE_CLOSE(0.32738, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the perfect cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyPerfect) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first perfect reactant (numI=5)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Perfect, 5);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Perfect == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Perfect
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (I1)
	auto secondReactant = (AlloyCluster *) network->get(Species::I, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-895996635.6, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -3.15972e+08, -3.30963e+08, -3.41479e+08,
			-3.49851e+08, -3.56921e+08, -3.15972e+08, -3.30963e+08,
			-3.41479e+08, -3.49851e+08, -3.63101e+08, -3.68627e+08,
			-3.73648e+08, -3.78267e+08, -3.82556e+08, -3.86567e+08,
			-3.90342e+08, -3.93912e+08, -3.97304e+08, -4.00538e+08,
			-4.03631e+08, -4.06598e+08, -4.09451e+08, -4.12199e+08,
			-4.14853e+08, -3.9251e+08, -4.03264e+08, -4.13273e+08, -4.22674e+08,
			-4.31565e+08, -4.40023e+08, -4.48103e+08, -4.55854e+08,
			-4.63311e+08, -4.70507e+08, -4.77466e+08, -4.84212e+08,
			-4.90761e+08, -4.97131e+08, -5.03336e+08, -3.26801e+09,
			-6.95882e+08, -6.67353e+08, -6.47292e+08, -6.32844e+08,
			-6.22294e+08, -6.14552e+08, -6.08892e+08, -6.04813e+08,
			-6.01958e+08, -6.00067e+08, 0, 0, 0, 0, 0, -3.80816e+08,
			-3.9251e+08, -4.03264e+08, -4.13273e+08, -4.22674e+08, -4.31565e+08,
			-4.40023e+08, -4.48103e+08, -4.55854e+08, -4.63311e+08,
			-4.70507e+08, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Perfect5
	BOOST_REQUIRE_CLOSE(0.270047, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the frank cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyFrank) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first frank reactant (numI=5)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Frank, 5);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Frank == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Frank
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Faulted
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			// Frank
			1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (I1)
	auto secondReactant = (AlloyCluster *) network->get(Species::I, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-0.007587566, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -4.50078e-06, -4.70662e-06, -4.85102e-06,
			-4.96597e-06, -5.06304e-06, -0.0151751, -0.0158692, -0.016356,
			0.0107396, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.80816e+08, -3.26152e+08,
			-2.865e+08, -2.5634e+08, -2.32577e+08, -2.13337e+08, -1.97416e+08,
			-1.84007e+08, -1.72546e+08, -1.62628e+08, -1.53954e+08, 0, 0, 0, 0,
			0, -0.0151751, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Frank5
	BOOST_REQUIRE_CLOSE(0.2988566, reactant->getReactionRadius(), 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
