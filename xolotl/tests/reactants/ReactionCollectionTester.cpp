#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <experimental/FeReactionNetwork.h>
#include <experimental/NEReactionNetwork.h>
#include <experimental/PSIReactionNetwork.h>

using namespace std;
using namespace xolotlCore;

class KokkosContext {
public:
	KokkosContext() {
		::Kokkos::initialize();
	}

	~KokkosContext() {
		::Kokkos::finalize();
	}
};
BOOST_GLOBAL_FIXTURE(KokkosContext);

BOOST_AUTO_TEST_SUITE(ReactionCollection_testSuite)

BOOST_AUTO_TEST_CASE(FeNetwork) {
    using NetworkType = experimental::FeReactionNetwork;
    using ReactionCollection =
        experimental::detail::ReactionCollection<NetworkType>;

    using ProdReaction = typename NetworkType::Traits::ProductionReactionType;
    using DissReaction = typename NetworkType::Traits::DissociationReactionType;
    using SinkReaction = typename NetworkType::Traits::SinkReactionType;
    using ReSoReaction = typename NetworkType::Traits::ReSolutionReactionType;

    Kokkos::View<ProdReaction*> prodReactions("prod", 3);
    Kokkos::View<DissReaction*> dissReactions("diss", 4);
    Kokkos::View<SinkReaction*> sinkReactions("sink", 5);
    Kokkos::View<ReSoReaction*> resoReactions("reso", 6);

    ReactionCollection rColl;
    rColl.setView(dissReactions);
    rColl.setView(prodReactions);
    rColl.setView(Kokkos::View<DissReaction*>{});
    // auto rColl = ReactionCollection(prodReactions, dissReactions, sinkReactions,
    //     resoReactions);
    std::cout << rColl.getDeviceMemorySize() << '\n';
    rColl.apply(DEVICE_LAMBDA (auto&& reaction) {
        std::cout << typeid(reaction).name() << '\n';
    });
    int numReactions = 0;
    rColl.reduce(DEVICE_LAMBDA (auto&& reaction, int& local) {
        ++local;
    }, numReactions);
    Kokkos::fence();
    // BOOST_REQUIRE(numReactions == prodReactions.size() + dissReactions.size() +
    //     sinkReactions.size() + resoReactions.size());

    // auto prodView = rColl.getView<ProdReaction>();
    // BOOST_REQUIRE(prodView.size() == prodReactions.size());
    // auto resoView = rColl.getView<ReSoReaction>();
    // BOOST_REQUIRE(resoView.size() == resoReactions.size());
}

BOOST_AUTO_TEST_SUITE_END()
