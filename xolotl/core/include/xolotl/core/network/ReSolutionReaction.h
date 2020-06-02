#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
class ReSolutionReaction : public Reaction<TNetwork, TDerived>
{
    friend class Reaction<TNetwork, TDerived>;

public:
    using NetworkType = TNetwork;
    using Superclass = Reaction<TNetwork, TDerived>;
    using ClusterDataRef = typename Superclass::ClusterDataRef;
    using IndexType = typename Superclass::IndexType;
    using Connectivity = typename Superclass::Connectivity;
    using ConcentrationsView = typename Superclass::ConcentrationsView;
    using FluxesView = typename Superclass::FluxesView;
    using AmountType = typename Superclass::AmountType;
    using ReactionDataRef = typename Superclass::ReactionDataRef;

    ReSolutionReaction() = default;

    KOKKOS_INLINE_FUNCTION
    ReSolutionReaction(ReactionDataRef reactionData,
        ClusterDataRef clusterData, IndexType reactionId,
        IndexType cluster0, IndexType cluster1, IndexType cluster2);

    KOKKOS_INLINE_FUNCTION
    ReSolutionReaction(ReactionDataRef reactionData,
        ClusterDataRef clusterData, IndexType reactionId,
        const detail::ClusterSet& clusterSet);

    static
    detail::CoefficientsView
    allocateCoefficientsView(IndexType size)
    {
        return detail::CoefficientsView("ReSolution Coefficients", size,
            Superclass::coeffsSingleExtent, 1, 3,
            Superclass::coeffsSingleExtent);
    }

private:
    KOKKOS_INLINE_FUNCTION
    void
    computeCoefficients();

    KOKKOS_INLINE_FUNCTION
    double
    computeRate(IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    computeConnectivity(const Connectivity& connectivity);

    KOKKOS_INLINE_FUNCTION
    void
    computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    computePartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values, Connectivity connectivity,
        IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    double
    computeLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
        IndexType gridIndex);

protected:
    IndexType _reactant;
    static constexpr auto invalidIndex = Superclass::invalidIndex;
    Kokkos::Array<IndexType, 2> _products {invalidIndex, invalidIndex};

    static constexpr auto nMomentIds = Superclass::nMomentIds;
    Kokkos::Array<IndexType, nMomentIds> _reactantMomentIds;
    Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _productMomentIds;
};
}
}
}

#include <xolotl/core/network/ReSolutionReaction.inl>
#include <xolotl/core/network/detail/ReSolutionReactionGenerator.h>
