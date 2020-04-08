#pragma once

#include <experimental/ReactionNetwork.h>
#include <experimental/FeReaction.h>
#include <experimental/FeTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
class FeReactionGenerator;
}

class FeReactionNetwork :
    public ReactionNetwork<FeReactionNetwork>
{
    friend class ReactionNetwork<FeReactionNetwork>;
    friend class detail::ReactionNetworkWorker<FeReactionNetwork>;

public:
    using Superclass = ReactionNetwork<FeReactionNetwork>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using AmountType = typename Superclass::AmountType;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;
    using ClusterDataRef = typename Superclass::ClusterDataRef;

    using Superclass::Superclass;

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return tungstenLatticeConstant;
        }
        return latticeParameter;
    }

    double
    checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return heliumRadius;
        }
        return impurityRadius;
    }

    detail::FeReactionGenerator
    getReactionGenerator() const noexcept;
};

namespace detail
{
class FeReactionGenerator : public
    ReactionGenerator<FeReactionNetwork,
        FeReactionGenerator>
{
public:
    using Network = FeReactionNetwork;
    using Subpaving = typename Network::Subpaving;
    using IndexType = typename Network::IndexType;

    using Superclass = ReactionGenerator<FeReactionNetwork,
        FeReactionGenerator>;

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
}
}
}

#include <experimental/FeClusterGenerator.h>
#include <experimental/FeReactionNetwork.inl>
