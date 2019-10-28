#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(AmountType maxSpeciesAmount)
    :
    _subpaving(Region{{
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1}}},
    {{{
        maxSpeciesAmount+1,
        maxSpeciesAmount+1,
        maxSpeciesAmount+1,
        maxSpeciesAmount+1,
        maxSpeciesAmount+1}}})
{
}


template <typename TImpl>
typename ReactionNetwork<TImpl>::Cluster
ReactionNetwork<TImpl>::getCluster(const Composition& comp)
{
    _subpaving.syncAll(plsm::onHost);
    Cluster ret(*this, _subpaving.findTileId(comp, plsm::onHost));
    return ret;
}


template <typename TImpl>
void
ReactionNetwork<TImpl>::defineMomentIds()
{
    constexpr auto invalid = plsm::invalid<std::size_t>;
    _subpaving.syncAll(plsm::onHost);
    auto tiles = _subpaving.getTiles(plsm::onHost);
    auto nClusters = _subpaving.getNumberOfTiles(plsm::onHost);

    //FIXME: _momentIds data lives on the device this way
    _momentIds = Kokkos::View<std::size_t*[4]>("Moment Ids", nClusters);

    auto current = nClusters;
    for (std::size_t c = 0; c < nClusters; ++c) {
        const auto& reg = tiles(c).getRegion();
        for (auto k : getSpeciesRangeNoI()) {
            if (reg[k].length() == 1) {
                _momentIds(c, k()) = invalid;
            }
            else {
                _momentIds(c, k()) = current;
                ++current;
            }
        }
    }
}
}
}
