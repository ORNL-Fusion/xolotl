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
template <typename TSpeciesEnum>
class FeReactionGenerator;
}

template <typename TSpeciesEnum>
class FeReactionNetwork :
    public ReactionNetwork<FeReactionNetwork<TSpeciesEnum>>
{
    friend class ReactionNetwork<FeReactionNetwork<TSpeciesEnum>>;
    friend class detail::ReactionNetworkWorker<FeReactionNetwork<TSpeciesEnum>>;

public:
    using Superclass = ReactionNetwork<FeReactionNetwork<TSpeciesEnum>>;
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

    void
    addReactionFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex)
    {
        return;
    }

    void
    addReactionPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, IndexType gridIndex)
    {
        return;
    }

    detail::FeReactionGenerator<Species>
    getReactionGenerator() const noexcept
    {
        return detail::FeReactionGenerator<Species>{*this};
    }
};

namespace detail
{
template <typename TSpeciesEnum>
class FeReactionGenerator : public
    ReactionGenerator<FeReactionNetwork<TSpeciesEnum>,
        FeReactionGenerator<TSpeciesEnum>>
{
public:
    using Network = FeReactionNetwork<TSpeciesEnum>;
    using Subpaving = typename Network::Subpaving;
    using IndexType = typename Network::IndexType;

    using Superclass = ReactionGenerator<FeReactionNetwork<TSpeciesEnum>,
        FeReactionGenerator<TSpeciesEnum>>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;
};
}
}
}

#include <experimental/FeClusterGenerator.h>
#include <experimental/FeReactionNetwork.inl>
