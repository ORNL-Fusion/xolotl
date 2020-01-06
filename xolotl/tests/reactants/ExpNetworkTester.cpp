#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <unordered_map>
#include <experimental/PSIReactionNetwork.h>

using namespace std;
using namespace xolotlCore;

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
BOOST_GLOBAL_FIXTURE(KokkosContext);

/**
 * This suite is responsible for testing the experimental networks.
 */
BOOST_AUTO_TEST_SUITE(ExpNetwork_testSuite)

/**
 * This operation checks the creation of a PSI network
 */
BOOST_AUTO_TEST_CASE(PSINetwork) {
	std::cout <<"in test" << std::endl;

    using NetworkType =
        experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
    using Spec = typename NetworkType::Species;
    using Composition = typename NetworkType::Composition;

    // auto rNetwork = experimental::makeSimpleReactionNetwork<NetworkType>();
    NetworkType rNetwork({8, 1, 1, 2, 6}, 0, Options{});
	std::cout <<"created network" << std::endl;

	rNetwork.syncClusterDataOnHost();
    Composition comp{};
    comp[Spec::He] = 3;
    comp[Spec::V] = 2;
    comp[Spec::I] = 0;
    auto cluster = rNetwork.findCluster(comp, plsm::onHost);
    auto compRegion = cluster.getRegion();
    BOOST_REQUIRE(compRegion.isSimplex());
    BOOST_REQUIRE_EQUAL(compRegion.getOrigin(), comp);

    cluster = rNetwork.getCluster(0, plsm::onHost);
    compRegion = cluster.getRegion();
    BOOST_REQUIRE_EQUAL(compRegion.getOrigin(), Composition{});

    using Reaction = typename NetworkType::ReactionType;
    // Reaction r0;
    // Reaction r(rNetwork, 0, Reaction::Type::dissociation, 0, 1, 2);
    // BOOST_REQUIRE(r.getType() == Reaction::Type::dissociation);

    auto concs = new double[3]{};
    auto fluxes = new double[3]{};
    using Vue = Reaction::ConcentrationsView;
    // r.contributeFlux(Vue(concs, 3), Vue(fluxes, 3), 0);
}

BOOST_AUTO_TEST_SUITE_END()
