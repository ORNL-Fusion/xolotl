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
    getRegion() const
    {
        //FIXME: explicitly using host space
        return _network->_subpaving.getTiles(plsm::onHost)(_id).getRegion();
    }

    decltype(auto)
    getMomentIds()
    {
        return Kokkos::subview(_network->_momentIds, _id, Kokkos::ALL);
    }

    double
    getReactionRadius()
    {
        return _network->_reactionRadius(_id);
    }

    double
    getFormationEnergy()
    {
        return _network->_formationEnergy(_id);
    }

    double
    getDiffusionCoefficient(std::size_t gridIndex)
    {
        return _network->_diffusionCoefficient(_id, gridIndex);
    }

private:
    ReactionNetwork<TImpl>* _network {nullptr};
    std::size_t _id {plsm::invalid<std::size_t>};
};
}
}
