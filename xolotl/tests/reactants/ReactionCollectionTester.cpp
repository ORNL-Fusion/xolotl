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

    using ReactionTypes = std::tuple<ProdReaction, DissReaction, SinkReaction>;

    Kokkos::View<ProdReaction*> prodReactions("prod", 3);
    Kokkos::View<DissReaction*> dissReactions("diss", 4);
    Kokkos::View<SinkReaction*> sinkReactions("sink", 5);

    auto rColl = ReactionCollection(0, prodReactions, dissReactions, sinkReactions);
//    rColl.apply(DEVICE_LAMBDA (auto&& reaction) {
//        //FIXME: Can't use typeid in device code. Use an array of three and
//        //specialized structs to convert between reaction types and indices
//        std::cout << typeid(reaction).name() << '\n';
//    });
    int numReactions = 0;
    rColl.reduce(DEVICE_LAMBDA (auto&& reaction, int& local) {
        ++local;
    }, numReactions);
    Kokkos::fence();
    BOOST_REQUIRE(numReactions == prodReactions.size() + dissReactions.size() +
        sinkReactions.size());

    auto prodView = rColl.getView<ProdReaction>();
    BOOST_REQUIRE(prodView.size() == prodReactions.size());
    auto sinkView = rColl.getView<SinkReaction>();
    BOOST_REQUIRE(sinkView.size() == sinkReactions.size());
}

BOOST_AUTO_TEST_SUITE_END()
