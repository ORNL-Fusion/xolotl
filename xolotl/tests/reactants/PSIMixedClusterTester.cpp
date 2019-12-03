#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <PSIHeCluster.h>
#include <PSIMixedCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <math.h>
#include <limits>

#include <experimental/PSIReactionNetwork.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

class KokkosContext
{
public:
    KokkosContext()
    {
        ::Kokkos::initialize();
    }

    ~KokkosContext()
    {
        ::Kokkos::finalize();
    }
};
static KokkosContext kokkosContext{};

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the PSIMixedCluster.
 */
BOOST_AUTO_TEST_SUITE(PSIMixedCluster_testSuite)

BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	// Create a simple reaction network and create a Mixed cluster
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork(0);
	PSIMixedCluster cluster(4, 0, 0, 5, *(network.get()), registry);

	// Get the composition back
	auto composition = cluster.getComposition();

	// Check the composition is the created one
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 4);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 5);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);

	// Check if it is a mixed cluster
	BOOST_REQUIRE_EQUAL(cluster.isMixed(), true);

	return;
}

/**
 * This operation checks the ability of the PSIMixedCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

    using NetworkType =
        experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
    using Spec = typename NetworkType::Species;
    using Composition = typename NetworkType::Composition;

    // auto rNetwork = experimental::makeSimpleReactionNetwork<NetworkType>();
    NetworkType rNetwork({10, 10, 10, 10, 10}, 0, Options{});
    Composition comp{};
    comp[Spec::He] = 3;
    comp[Spec::V] = 2;
    comp[Spec::I] = 0;
    auto cluster = rNetwork.findCluster(comp);
    auto compRegion = cluster.getRegion();
    BOOST_REQUIRE(compRegion.isSimplex());
    BOOST_REQUIRE_EQUAL(compRegion.getOrigin(), comp);

    cluster = rNetwork.getCluster(0);
    compRegion = cluster.getRegion();
    BOOST_REQUIRE_EQUAL(compRegion.getOrigin(), Composition{});

    using Reaction = typename NetworkType::ReactionType;
    Reaction r0;
    Reaction r(rNetwork, 0, Reaction::Type::dissociation, 0, 1, 2);
    BOOST_REQUIRE(r.getType() == Reaction::Type::dissociation);

    auto concs = new double[3]{};
    auto fluxes = new double[3]{};
    using Vue = Reaction::ConcentrationsView;
    // r.contributeFlux(Vue(concs, 3), Vue(fluxes, 3), 0);

	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork();

	// Check the reaction connectivity of the PSIMixed cluster
	// with 3He and 2V
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 3;
	composition[toCompIdx(Species::V)] = 2;
	composition[toCompIdx(Species::I)] = 0;
	auto reactant = (PSICluster *) network->get(ReactantType::PSIMixed,
			composition);

	// Check the type name
	BOOST_REQUIRE(ReactantType::PSIMixed == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 3);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 2);

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0,

			// D
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// T
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// V
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0,

			// I
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0,

			// Mixed
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1,
			1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// HeI
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the ability of the PSIMixedCluster to compute the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {
	// Local Declarations
	auto network = getSimplePSIReactionNetwork();
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an Mixed cluster with compostion 2,1,0.
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 2;
	composition[toCompIdx(Species::V)] = 1;
	composition[toCompIdx(Species::I)] = 0;
	auto cluster = (PSICluster *) network->get(ReactantType::PSIMixed,
			composition);
	// Get one that it combines with (He)
	auto secondCluster = (PSICluster *) network->get(Species::He, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the preprocessor.
	cluster->setDiffusionFactor(0.0);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(2.950E+10);
	secondCluster->setMigrationEnergy(0.13);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-1134677704810.4, flux, 0.1);

	return;
}

/**
 * This operation checks the PSIMixedCluster get*PartialDerivatives methods.
 */
// BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
// 	// Local Declarations
// 	// The vector of partial derivatives to compare with
// 	double knownPartials[] = { 0, 0, -2.74742, 0, -2.74742, 0, 0, 0, -1.85429,
// 			0, -689.98, 344.99, 344.99, 0, 0, 0 };
// 	// Get the simple reaction network
// 	auto network = getSimplePSIReactionNetwork(2);
// 	// Add a grid point for the rates
// 	network->addGridPoints(1);

// 	// Get an Mixed cluster with compostion 1,1,0.
// 	IReactant::Composition composition;
// 	composition[toCompIdx(Species::He)] = 1;
// 	composition[toCompIdx(Species::V)] = 1;
// 	composition[toCompIdx(Species::I)] = 0;
// 	auto cluster = (PSICluster *) network->get(ReactantType::PSIMixed,
// 			composition);
// 	// Set the diffusion factor and migration energy based on the
// 	// values from the tungsten benchmark for this problem.
// 	cluster->setDiffusionFactor(0.0);
// 	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
// 	cluster->setConcentration(0.5);

// 	// Compute the rate constants that are needed for the partials
// 	network->setTemperature(1000.0, 0);
// 	// Get the vector of partial derivatives
// 	auto partials = cluster->getPartialDerivatives(0);

// 	// Check the size of the partials
// 	BOOST_REQUIRE_EQUAL(partials.size(), 16U);

// 	// Check all the values
// 	for (unsigned int i = 0; i < partials.size(); i++) {
// 		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
// 	}

// 	return;
// }

/**
 * This operation checks the reaction radius for PSIMixedCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the Mixed cluster
	shared_ptr<PSIMixedCluster> cluster;

	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1372650265, 0.1778340462, 0.2062922619,
			0.2289478080, 0.2480795532 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<PSIMixedCluster>(
				new PSIMixedCluster(1, 0, 0, i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.000001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
