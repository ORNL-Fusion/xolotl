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
    using ReactionType = typename Superclass::ReactionType;
    using Composition = typename Superclass::Composition;
    using Species = typename Superclass::Species;
    using AmountType = typename Superclass::AmountType;

    static constexpr auto invalid = plsm::invalid<std::size_t>;

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
public:
    using Network = PSIReactionNetwork<TSpeciesEnum>;
    using Subpaving = typename Network::Subpaving;
    using ClusterSet = typename Network::ReactionType::ClusterSet;

    using Superclass = ReactionGenerator<PSIReactionNetwork<TSpeciesEnum>,
        PSIReactionGenerator<TSpeciesEnum>>;
    using Superclass::Superclass;

    template <typename TTag>
    KOKKOS_INLINE_FUNCTION
    void
    operator()(std::size_t i, std::size_t j, TTag tag) const;
};
}
}
}

#include <experimental/PSIClusterGenerator.h>
#include <experimental/PSIReactionNetwork.inl>
