#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
class ReactionNetwork<TImpl>::Cluster
{
public:
    using NetworkType = ReactionNetwork<TImpl>;
    using Region = typename NetworkType::Region;

    Cluster() = delete;
    Cluster(ReactionNetwork& network, std::size_t id)
        :
        _network{&network},
        _id{id}
    {
    }

    Region
    getRegion()
    {
        //FIXME: explicitly using host space
        return _network->getSubpaving().getTiles(plsm::onHost)(_id).getRegion();
    }

    double
    getReactionRadius()
    {
        return _network->getReactionRadius(_id);
    }

private:
    ReactionNetwork<TImpl>* _network {nullptr};
    std::size_t _id {plsm::invalid<std::size_t>};
};
}
}
