#pragma once

#include <xolotl/core/reactants/ReactionNetwork.h>
#include <xolotl/core/reactants/AlloyReaction.h>
#include <xolotl/core/reactants/AlloyTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
class AlloyReactionGenerator;

class AlloyClusterUpdater;
}

class AlloyReactionNetwork : public ReactionNetwork<AlloyReactionNetwork>
{
    friend class ReactionNetwork<AlloyReactionNetwork>;
    friend class detail::ReactionNetworkWorker<AlloyReactionNetwork>;

public:
    using Superclass = ReactionNetwork<AlloyReactionNetwork>;
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
            return alloyLatticeConstant;
        }
        return latticeParameter;
    }

    double checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return alloyCoreRadius;
        }
        return impurityRadius;
    }

    void
    checkTiles(const IOptions& options)
    {
        return;
    }

    detail::AlloyReactionGenerator
    getReactionGenerator() const noexcept;
};

namespace detail
{
class AlloyReactionGenerator :
    public ReactionGenerator<AlloyReactionNetwork, AlloyReactionGenerator>
{
public:
    using Network = AlloyReactionNetwork;
    using Subpaving = typename Network::Subpaving;
    using Superclass =
        ReactionGenerator<AlloyReactionNetwork, AlloyReactionGenerator>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    addSinks(IndexType i, TTag tag) const;

    ReactionCollection<Network>
    getReactionCollection() const;
};

class AlloyClusterUpdater
{
public:
    using Network = AlloyReactionNetwork;
    using ClusterData = typename Network::ClusterData;
    using IndexType = typename Network::IndexType;
};
}
}
}

#include <xolotl/core/reactants/AlloyClusterGenerator.h>
#include <xolotl/core/reactants/AlloyReactionNetwork.inl>
