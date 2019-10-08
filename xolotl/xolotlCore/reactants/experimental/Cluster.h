#pragma once

#include <array>
#pragma once

namespace xolotlCore
{
namespace experimental
{
template <std::size_t NumSpecies>
class ReactionNetwork<NumSpecies>::Cluster
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
        return _network->getSubpaving().getTilesOnHost()(_id).getRegion();
    }

private:
    ReactionNetwork<NumSpecies>* _network {nullptr};
    std::size_t _id {plsm::invalid<std::size_t>};
};
}
}
