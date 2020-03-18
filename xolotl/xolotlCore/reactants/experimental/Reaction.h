#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/detail/KokkosExtension.h>
#include <plsm/Utility.h>

#include <experimental/Cluster.h>
#include <experimental/IReactionNetwork.h>
#include <experimental/ReactionData.h>
#include <experimental/ReactionNetworkTraits.h>
#include <experimental/SpeciesEnumSequence.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
struct ClusterSet
{
    using IndexType = ReactionNetworkIndexType;
    static constexpr IndexType invalidIndex = IReactionNetwork::invalidIndex();

    IndexType cluster0 {invalidIndex};
    IndexType cluster1 {invalidIndex};
    IndexType cluster2 {invalidIndex};
    IndexType cluster3 {invalidIndex};

    ClusterSet() = default;

    KOKKOS_INLINE_FUNCTION
    ClusterSet(IndexType cl0, IndexType cl1, IndexType cl2 = invalidIndex,
            IndexType cl3 = invalidIndex)
        :
        cluster0{cl0},
        cluster1{cl1},
        cluster2{cl2},
        cluster3{cl3}
    {
    }

    KOKKOS_INLINE_FUNCTION
    bool
    valid() const noexcept
    {
        return cluster0 != invalidIndex;
    }
};
}

template <typename TNetwork, typename TDerived>
class Reaction
{
    using Types = detail::ReactionNetworkTypes<TNetwork>;
    using Props = detail::ReactionNetworkProperties<TNetwork>;

protected:
    static constexpr auto invalidIndex = detail::InvalidIndex::value;
    static constexpr auto nMomentIds = Props::numSpeciesNoI;

public:
    using NetworkType = TNetwork;
    using Species = typename Types::Species;
    using IndexType = typename Types::IndexType;
    using AmountType = typename Types::AmountType;
    using Region = typename Types::Region;
    using Composition = typename Types::Composition;
    using ClusterDataRef = typename Types::ClusterDataRef;
    using ConcentrationsView = IReactionNetwork::ConcentrationsView;
    using FluxesView = IReactionNetwork::FluxesView;
    using Connectivity = typename IReactionNetwork::Connectivity;

    Reaction() = default;

    KOKKOS_INLINE_FUNCTION
    Reaction(detail::ReactionDataRef reactionData, ClusterDataRef clusterData,
        IndexType reactionId);

    KOKKOS_INLINE_FUNCTION
    void
    updateData(detail::ReactionDataRef reactionData, ClusterDataRef clusterData,
        IndexType reactionId);

    KOKKOS_INLINE_FUNCTION
    void
    updateRates()
    {
        for (IndexType i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeRate(i);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void
    contributeConnectivity(const Connectivity& connectivity)
    {
        asDerived()->computeConnectivity(connectivity);
    }

    KOKKOS_INLINE_FUNCTION
    void
    contributeFlux(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex)
    {
        asDerived()->computeFlux(concentrations, fluxes, gridIndex);
    }

    KOKKOS_INLINE_FUNCTION
    void
    contributePartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values, Connectivity connectivity,
        IndexType gridIndex)
    {
        asDerived()->computePartialDerivatives(concentrations, values,
            connectivity, gridIndex);
    }

    KOKKOS_INLINE_FUNCTION
    double
    contributeLeftSideRate(ConcentrationsView concentrations,
        IndexType clusterId, IndexType gridIndex)
    {
        return asDerived()->computeLeftSideRate(concentrations, clusterId,
            gridIndex);
    }

protected:
    KOKKOS_INLINE_FUNCTION
    TDerived*
    asDerived()
    {
        return static_cast<TDerived*>(this);
    }

    KOKKOS_INLINE_FUNCTION
    void
    initialize()
    {
        asDerived()->computeCoefficients();
        updateRates();
    }

    KOKKOS_INLINE_FUNCTION
    AmountType
    computeOverlap(const Region& singleClReg, const Region& pairCl1Reg,
        const Region& pairCl2Reg);

    KOKKOS_INLINE_FUNCTION
    void
    copyMomentIds(IndexType clusterId,
        Kokkos::Array<IndexType, nMomentIds>& momentIds)
    {
        if (clusterId == invalidIndex) {
            for (IndexType i = 0; i < nMomentIds; ++i) {
                momentIds[i] = invalidIndex;
            }
            return;
        }

        const auto& mIds = _clusterData.getCluster(clusterId).getMomentIds();
        for (IndexType i = 0; i < nMomentIds; ++i) {
            momentIds[i] = mIds[i];
        }
    }

    KOKKOS_INLINE_FUNCTION
    void
    addConnectivity(IndexType rowId, IndexType columnId,
        const Connectivity& connectivity)
    {
        connectivity.add(rowId, columnId);
    }

protected:
    ClusterDataRef _clusterData;

    //! Reaction rate (k)
    using RateSubView = decltype(
        std::declval<detail::ReactionDataRef>().getRates(0));
    RateSubView _rate;

    //! Flux coefficients
    using CoefsSubView = decltype(
        std::declval<detail::ReactionDataRef>().getCoefficients(0));
    CoefsSubView _coefs;
};

template <typename TNetwork, typename TDerived>
class ProductionReaction : public Reaction<TNetwork, TDerived>
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
    using Composition = typename Superclass::Composition;
    using Region = typename Superclass::Region;
    using AmountType = typename Superclass::AmountType;

    ProductionReaction() = default;

    KOKKOS_INLINE_FUNCTION
    ProductionReaction(detail::ReactionDataRef reactionData,
        ClusterDataRef clusterData, IndexType reactionId,
        IndexType cluster0, IndexType cluster1,
        IndexType cluster2 = invalidIndex, IndexType cluster3 = invalidIndex);

    KOKKOS_INLINE_FUNCTION
    ProductionReaction(detail::ReactionDataRef reactionData,
        ClusterDataRef clusterData, IndexType reactionId,
        const detail::ClusterSet& clusterSet);

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
    static constexpr auto invalidIndex = Superclass::invalidIndex;
    Kokkos::Array<IndexType, 2> _reactants {invalidIndex, invalidIndex};
    Kokkos::Array<IndexType, 2> _products {invalidIndex, invalidIndex};

    static constexpr auto nMomentIds = Superclass::nMomentIds;
    Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _reactantMomentIds;
    Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _productMomentIds;
};

template <typename TNetwork, typename TDerived>
class DissociationReaction : public Reaction<TNetwork, TDerived>
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

    DissociationReaction() = default;

    KOKKOS_INLINE_FUNCTION
    DissociationReaction(detail::ReactionDataRef reactionData,
        ClusterDataRef clusterData, IndexType reactionId,
        IndexType cluster0, IndexType cluster1, IndexType cluster2);

    KOKKOS_INLINE_FUNCTION
    DissociationReaction(detail::ReactionDataRef reactionData,
        ClusterDataRef clusterData, IndexType reactionId,
        const detail::ClusterSet& clusterSet);

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

#include <experimental/Reaction.inl>
