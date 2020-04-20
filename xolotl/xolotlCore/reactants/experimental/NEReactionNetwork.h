#pragma once

#include <experimental/ReactionNetwork.h>
#include <experimental/NEReaction.h>
#include <experimental/NETraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
class NEReactionGenerator;

class NEClusterUpdater;
}

class NEReactionNetwork : public ReactionNetwork<NEReactionNetwork>
{
    friend class ReactionNetwork<NEReactionNetwork>;
    friend class detail::ReactionNetworkWorker<NEReactionNetwork>;

public:
    using Superclass = ReactionNetwork<NEReactionNetwork>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;

    using Superclass::Superclass;

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return uraniumDioxydeLatticeConstant;
        }
        return latticeParameter;
    }

    double checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return xenonRadius;
        }
        return impurityRadius;
    }

    detail::NEReactionGenerator
    getReactionGenerator() const noexcept;
};

namespace detail
{
class NEReactionGenerator :
    public ReactionGenerator<NEReactionNetwork, NEReactionGenerator>
{
public:
    using Network = NEReactionNetwork;
    using Subpaving = typename Network::Subpaving;
    using Superclass =
        ReactionGenerator<NEReactionNetwork, NEReactionGenerator>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    addSinks(IndexType i, TTag tag) const;
};

class NEClusterUpdater
{
public:
    using Network = NEReactionNetwork;
    using ClusterData = typename Network::ClusterData;
    using IndexType = typename Network::IndexType;

    KOKKOS_INLINE_FUNCTION
    void
    updateDiffusionCoefficient(const ClusterData& data, IndexType clusterId,
        IndexType gridIndex) const
    {
        data.diffusionCoefficient(clusterId, gridIndex) =
            data.diffusionFactor(clusterId) * exp(
                -data.migrationEnergy(clusterId) /
                (kBoltzmann * data.temperature(gridIndex)));
    }
};
}
}
}

#include <experimental/NEClusterGenerator.h>
#include <experimental/NEReactionNetwork.inl>
