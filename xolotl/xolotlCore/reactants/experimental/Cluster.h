#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
class ReactionNetwork<TImpl>::Cluster
{
public:
    Cluster() = delete;
    Cluster(ReactionNetwork& network, std::size_t id)
        :
        _network{&network},
        _id{id}
    {
    }

    decltype(auto)
    getRegion()
    {
        return _network->getSubpaving().getTiles(plsm::onHost)(_id).getRegion();
    }

private:
    ReactionNetwork<TImpl>* _network {nullptr};
    std::size_t _id {plsm::invalid<std::size_t>};
};
}
}
