#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/detail/KokkosExtension.h>
#include <plsm/Utility.h>

#include <experimental/Cluster.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
struct ReactionDataRef
{
    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getCoefficients(std::size_t reactionId)
    {
        if (reactionId < productionCoeffs.extent(0)) {
            return Kokkos::subview(productionCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
        else {
            reactionId -= productionCoeffs.extent(0);
            return Kokkos::subview(dissociationCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
    }

    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getRates(std::size_t reactionId)
    {
        return Kokkos::subview(rates, reactionId, Kokkos::ALL);
    }

    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> productionCoeffs;
    Kokkos::View<double*****, Kokkos::MemoryUnmanaged> dissociationCoeffs;
    Kokkos::View<double**, Kokkos::MemoryUnmanaged> rates;

    Kokkos::View<std::size_t**, Kokkos::MemoryUnmanaged> inverseMap;
};
}

template <typename TImpl>
struct ReactionNetworkTraits
{
};

namespace detail
{
template <typename TImpl>
struct ReactionNetworkTypes
{
    using AmountType = std::uint32_t;
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;
    using Subpaving = plsm::Subpaving<AmountType, Traits::numSpecies, Species>;
    using Region = typename Subpaving::RegionType;
    using Composition = typename Subpaving::PointType;
    using ClusterData = detail::ClusterData<Subpaving>;
    using ClusterDataMirror = detail::ClusterData<Subpaving, plsm::OnHost>;
    using ClusterDataRef = detail::ClusterDataRef<Subpaving>;
    using ConnectivityView = Kokkos::View<size_t**, Kokkos::MemoryUnmanaged>;
};

template <typename TImpl>
struct ReactionNetworkProperties
{
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;
    static constexpr std::size_t numSpecies = Traits::numSpecies;
    using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
    static constexpr std::size_t numSpeciesNoI = SpeciesSequence::sizeNoI();
};
}

template <typename TNetwork, typename TDerived>
class Reaction
{
    static constexpr std::size_t invalid = plsm::invalid<std::size_t>;

    using Types = detail::ReactionNetworkTypes<TNetwork>;
    using Props = detail::ReactionNetworkProperties<TNetwork>;

public:
    using NetworkType = TNetwork;
    using AmountType = typename Types::AmountType;
    using Region = typename Types::Region;
    using Composition = typename Types::Composition;
    using ClusterDataRef = typename Types::ClusterDataRef;
    using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using ConnectivityView = typename Types::ConnectivityView;

    enum class Type
    {
        production,
        dissociation
    };

    struct ClusterSet
    {
        std::size_t cluster0 {invalid};
        std::size_t cluster1 {invalid};
        std::size_t cluster2 {invalid};
        std::size_t cluster3 {invalid};

        ClusterSet() = default;

        ClusterSet(std::size_t cl0, std::size_t cl1, std::size_t cl2 = invalid,
                std::size_t cl3 = invalid)
            :
            cluster0{cl0},
            cluster1{cl1},
            cluster2{cl2},
            cluster3{cl3}
        {
        }
    };

    Reaction() = default;

    KOKKOS_INLINE_FUNCTION
    Reaction(detail::ReactionDataRef reactionData, ClusterDataRef clusterData,
        std::size_t reactionId, Type reactionType,
        std::size_t cluster0, std::size_t cluster1,
        std::size_t cluster2 = invalid, std::size_t cluster3 = invalid);

    KOKKOS_INLINE_FUNCTION
    Type
    getType() const noexcept
    {
        return _type;
    }

    KOKKOS_INLINE_FUNCTION
    void
    updateRates();

    KOKKOS_INLINE_FUNCTION
    void
    productionConnectivity(ConnectivityView connectivity);

    KOKKOS_INLINE_FUNCTION
    void
    dissociationConnectivity(ConnectivityView connectivity);

    KOKKOS_INLINE_FUNCTION
    void
    contributeConnectivity(ConnectivityView connectivity)
    {
        _connectFn(*this, connectivity);
    }

    KOKKOS_INLINE_FUNCTION
    void
    productionFlux(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    dissociationFlux(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    contributeFlux(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex)
    {
        _fluxFn(*this, concentrations, fluxes, gridIndex);
    }

    KOKKOS_INLINE_FUNCTION
    void
    productionPartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values,
        std::size_t gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    dissociationPartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values,
        std::size_t gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    contributePartialDerivatives(ConcentrationsView concentrations,
        Kokkos::View<double*> values,
        std::size_t gridIndex)
    {
        _partialsFn(*this, concentrations, values, gridIndex);
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
    computeOverlap(const Region& singleCl, const Region& pairCl1,
        const Region& pairCl2);

    KOKKOS_INLINE_FUNCTION
    void
    computeProductionCoefficients();

    KOKKOS_INLINE_FUNCTION
    void
    computeDissociationCoefficients();

    KOKKOS_INLINE_FUNCTION
    void
    copyMomentIds(std::size_t clusterId,
        Kokkos::Array<std::size_t, 4>& momentIds)
    {
        if (clusterId == invalid) {
            momentIds = {invalid, invalid, invalid, invalid};
            return;
        }

        const auto& mIds = _clusterData.getCluster(clusterId).getMomentIds();
        for (std::size_t i = 0; i < 4; ++i) {
            momentIds[i] = mIds[i];
        }
    }

    KOKKOS_INLINE_FUNCTION
    double
    computeProductionRate(std::size_t gridIndex);

    KOKKOS_INLINE_FUNCTION
    double
    computeDissociationRate(std::size_t gridIndex);

    KOKKOS_INLINE_FUNCTION
    void
    computeProductionRates()
    {
        for (std::size_t i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeProductionRate(i);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void
    computeDissociationRates()
    {
        for (std::size_t i = 0; i < _rate.extent(0); ++i) {
            _rate(i) = asDerived()->computeDissociationRate(i);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void
    addConnectivity(std::size_t rowId, std::size_t columnId,
        ConnectivityView connectivity)
    {
        // Check that the Ids are valid
        if (rowId == invalid || columnId == invalid) return;

        // Check if the columnId is already present
        std::size_t i = 0;
        for (i = 0; i < connectivity.extent(1); ++i) {
            if (connectivity(rowId, i) == columnId) {
                // It is already present, don't do anything
                return;
            }
            if (connectivity(rowId, i) == invalid) {
                // We reached the end of the list and it was not present
                break;
            }
        }
        // Add the value
        connectivity(rowId, i) = columnId;
    }

protected:
    ClusterDataRef _clusterData;

    Type _type {};

    //TODO: This might not compile for device code
    using ConnectFn = decltype(std::mem_fn(&Reaction::contributeConnectivity));
    ConnectFn _connectFn {nullptr};

    using FluxFn = decltype(std::mem_fn(&Reaction::contributeFlux));
    FluxFn _fluxFn {nullptr};

    using PartialsFn =
        decltype(std::mem_fn(&Reaction::contributePartialDerivatives));
    PartialsFn _partialsFn {nullptr};

    //Cluster indices for LHS and RHS of the reaction
    //Dissociation reactions always have 1 input and 2 outputs
    //Production reactions always have 2 inputs, but may have 0, 1, or 2 outputs
    Kokkos::Array<std::size_t, 2> _reactants {invalid, invalid};
    Kokkos::Array<std::size_t, 2> _products {invalid, invalid};

    static constexpr auto nMomentIds = Props::numSpeciesNoI;
    Kokkos::Array<Kokkos::Array<std::size_t, nMomentIds>, 2> _reactantMomentIds;
    Kokkos::Array<Kokkos::Array<std::size_t, nMomentIds>, 2> _productMomentIds;

    //! Reaction rate (k)
    Kokkos::View<double*, Kokkos::MemoryUnmanaged> _rate;

    //! Flux coefficients
    Kokkos::View<double****, Kokkos::MemoryUnmanaged> _coefs;

    Kokkos::View<std::size_t**, Kokkos::MemoryUnmanaged> _inverseMap;
};
}
}

#include <experimental/Reaction.inl>
