#pragma once

#include <experimental/ReactionNetwork.h>
#include <experimental/PSIReaction.h>
#include <experimental/PSITraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TSpeciesEnum>
class PSIReactionGenerator;
}

template <typename TSpeciesEnum>
class PSIReactionNetwork :
    public ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>
{
    friend class ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
    friend class detail::ReactionNetworkWorker<PSIReactionNetwork<TSpeciesEnum>>;

public:
    using Superclass = ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
    using Subpaving = typename Superclass::Subpaving;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using AmountType = typename Superclass::AmountType;
    using IndexType = typename Superclass::IndexType;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;

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
    checkTiles(const IOptions& options)
    {
        return;
    }

    detail::PSIReactionGenerator<Species>
    getReactionGenerator() const noexcept
    {
        return detail::PSIReactionGenerator<Species>{*this};
    }
};

namespace detail
{
template <typename TSpeciesEnum>
class PSIReactionGenerator : public
    ReactionGenerator<PSIReactionNetwork<TSpeciesEnum>,
        PSIReactionGenerator<TSpeciesEnum>>
{
    friend class ReactionGeneratorBase<PSIReactionNetwork<TSpeciesEnum>,
        PSIReactionGenerator<TSpeciesEnum>>;

public:
    using NetworkType = PSIReactionNetwork<TSpeciesEnum>;
    using Subpaving = typename NetworkType::Subpaving;
    using IndexType = typename NetworkType::IndexType;

    using Superclass = ReactionGenerator<PSIReactionNetwork<TSpeciesEnum>,
        PSIReactionGenerator<TSpeciesEnum>>;

    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(IndexType i, IndexType j, TTag tag) const;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    addSinks(IndexType i, TTag tag) const;

private:
    ReactionCollection<NetworkType>
    getReactionCollection() const;
};
}
}
}

#include <experimental/PSIClusterGenerator.h>
#include <experimental/PSIReactionNetwork.inl>
