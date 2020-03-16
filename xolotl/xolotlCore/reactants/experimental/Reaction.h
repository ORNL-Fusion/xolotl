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
template <typename TNetwork, typename TDerived>
class Reaction
{
    static constexpr auto invalidIndex = detail::InvalidIndex::value;

    using Types = detail::ReactionNetworkTypes<TNetwork>;
    using Props = detail::ReactionNetworkProperties<TNetwork>;

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

    enum class Type
    {
        production,
        dissociation
    };

    struct ClusterSet
    {
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

    Reaction() = default;

    KOKKOS_INLINE_FUNCTION
    Reaction(detail::ReactionDataRef reactionData, ClusterDataRef clusterData,
        IndexType reactionId, Type reactionType,
        IndexType cluster0, IndexType cluster1,
        IndexType cluster2 = invalidIndex, IndexType cluster3 = invalidIndex);

    KOKKOS_INLINE_FUNCTION
    Reaction(detail::ReactionDataRef reactionData, ClusterDataRef clusterData,
        IndexType reactionId, Type reactionType,
        const ClusterSet& clusterSet);

    KOKKOS_INLINE_FUNCTION
    Type
    getType() const noexcept
    {
        return _type;
    }

    KOKKOS_INLINE_FUNCTION
    void
    updateData(detail::ReactionDataRef reactionData, ClusterDataRef clusterData,
        IndexType reactionId);

    KOKKOS_INLINE_FUNCTION
    void
    updateRates();

    KOKKOS_INLINE_FUNCTION
    void
    productionConnectivity(const Connectivity& connectivity);

    KOKKOS_INLINE_FUNCTION
    void
    dissociationConnectivity(const Connectivity& connectivity);

    KOKKOS_INLINE_FUNCTION
    void
    contributeConnectivity(const Connectivity& connectivity)
    {
        ((*this).*(_connectFn))(connectivity);
    }

    KOKKOS_INLINE_FUNCTION
    void
    productionFlux(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    dissociationFlux(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    contributeFlux(ConcentrationsView concentrations, FluxesView fluxes,
        IndexType gridIndex)
    {
        ((*this).*(_fluxFn))(concentrations, fluxes, gridIndex);
    }

    KOKKOS_INLINE_FUNCTION
    void
    productionPartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values, Connectivity connectivity,
        IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    dissociationPartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values, Connectivity connectivity,
        IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    contributePartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values, Connectivity connectivity,
        IndexType gridIndex)
    {
        ((*this).*(_partialsFn))(concentrations, values, connectivity,
            gridIndex);
    }

    KOKKOS_INLINE_FUNCTION
    double
    productionLeftSideRate(ConcentrationsView concentrations,
        IndexType clusterId, IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    double
    dissociationLeftSideRate(ConcentrationsView concentrations,
        IndexType clusterId, IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    double
    contributeLeftSideRate(ConcentrationsView concentrations,
        IndexType clusterId, IndexType gridIndex)
    {
        return ((*this).*(_leftSideFn))(concentrations, clusterId, gridIndex);
    }

private:
    KOKKOS_INLINE_FUNCTION
    TDerived*
    asDerived()
    {
        return static_cast<TDerived*>(this);
    }

    KOKKOS_INLINE_FUNCTION
    AmountType
    computeOverlap(const Region& singleClReg, const Region& pairCl1Reg,
        const Region& pairCl2Reg);

    KOKKOS_INLINE_FUNCTION
    void
    computeProductionCoefficients();

    KOKKOS_INLINE_FUNCTION
    void
    computeDissociationCoefficients();

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
    double
    computeProductionRate(IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    double
    computeDissociationRate(IndexType gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    computeProductionRates()
    {
        for (IndexType i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeProductionRate(i);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void
    computeDissociationRates()
    {
        for (IndexType i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeDissociationRate(i);
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

    Type _type {};

    using ConnectFn =
        void (Reaction::*)(const Connectivity&);
    ConnectFn _connectFn {nullptr};

    using FluxFn =
        void (Reaction::*)(ConcentrationsView, FluxesView, IndexType);
    FluxFn _fluxFn {nullptr};

    using PartialsFn =
        void (Reaction::*)(ConcentrationsView, Kokkos::View<double*>,
            Connectivity, IndexType);
    PartialsFn _partialsFn {nullptr};

    using LeftSideFn =
        double (Reaction::*)(ConcentrationsView, IndexType, IndexType);
    LeftSideFn _leftSideFn {nullptr};

    //Cluster indices for LHS and RHS of the reaction
    //Dissociation reactions always have 1 input and 2 outputs
    //Production reactions always have 2 inputs, but may have 0, 1, or 2 outputs
    Kokkos::Array<IndexType, 2> _reactants {invalidIndex, invalidIndex};
    Kokkos::Array<IndexType, 2> _products {invalidIndex, invalidIndex};

    Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _reactantMomentIds;
    Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _productMomentIds;

    //! Reaction rate (k)
    using RateSubView = decltype(
        std::declval<detail::ReactionDataRef>().getRates(0));
    RateSubView _rate;

    //! Flux coefficients
    using CoefsSubView = decltype(
        std::declval<detail::ReactionDataRef>().getCoefficients(0));
    CoefsSubView _coefs;
};
}
}

#include <experimental/Reaction.inl>
