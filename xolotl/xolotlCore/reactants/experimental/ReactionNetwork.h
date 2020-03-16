#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Crs.hpp>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <IOptions.h>
#include <Options.h>

#include <experimental/Cluster.h>
#include <experimental/IReactionNetwork.h>
#include <experimental/Reaction.h>
#include <experimental/SpeciesEnumSequence.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TImpl>
struct ReactionNetworkWorker;

template <typename TImpl, typename TDerived>
class ReactionGenerator;
}

template <typename TImpl>
class ReactionNetwork : public IReactionNetwork
{
    friend class detail::ReactionNetworkWorker<TImpl>;
    template <typename, typename>
    friend class detail::ReactionGenerator;

public:
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;

private:
    using Types = detail::ReactionNetworkTypes<TImpl>;

    static constexpr std::size_t numSpecies = Traits::numSpecies;

public:
    using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
    using SpeciesRange = EnumSequenceRange<Species, numSpecies>;
    using ReactionType = typename Traits::ReactionType;
    using ClusterGenerator = typename Traits::ClusterGenerator;
    using AmountType = typename IReactionNetwork::AmountType;
    using Subpaving = typename Types::Subpaving;
    using SubdivisionRatio = plsm::SubdivisionRatio<numSpecies>;
    using Composition = typename Types::Composition;
    using Region = typename Types::Region;
    using Ival = typename Region::IntervalType;
    using ConcentrationsView = typename IReactionNetwork::ConcentrationsView;
    using FluxesView = typename IReactionNetwork::FluxesView;
    using SparseFillMap = typename IReactionNetwork::SparseFillMap;
    using ClusterData = typename Types::ClusterData;
    using ClusterDataMirror = typename Types::ClusterDataMirror;
    using ClusterDataRef = typename Types::ClusterDataRef;

    template <typename PlsmContext>
    using Cluster = Cluster<TImpl, PlsmContext>;

    ReactionNetwork() = default;

    ReactionNetwork(const Subpaving& subpaving, std::size_t gridSize,
        const IOptions& options);

    ReactionNetwork(const Subpaving& subpaving, std::size_t gridSize);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        std::size_t gridSize, const IOptions& options);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        std::size_t gridSize, const IOptions& options);

    KOKKOS_INLINE_FUNCTION
    static
    constexpr std::size_t
    getNumberOfSpecies() noexcept
    {
        return SpeciesSequence::size();
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr std::size_t
    getNumberOfSpeciesNoI() noexcept
    {
        return SpeciesSequence::sizeNoI();
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr SpeciesRange
    getSpeciesRange() noexcept
    {
        return SpeciesRange{};
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr SpeciesRange
    getSpeciesRangeNoI() noexcept
    {
        return SpeciesRange(SpeciesSequence::first(),
            SpeciesSequence::lastNoI());
    }

    void
    setLatticeParameter(double latticeParameter) override;

    void
    setImpurityRadius(double impurityRadius) noexcept override
    {
        this->_impurityRadius =
            asDerived()->checkImpurityRadius(impurityRadius);
    }

    void
    setGridSize(std::size_t gridSize) override;

    void
    setTemperatures(const std::vector<double>& gridTemperatures) override;

    std::uint64_t
    getDeviceMemorySize() const noexcept override
    {
        std::uint64_t ret = _subpaving.getDeviceMemorySize();

        ret += sizeof(_clusterData.numClusters);
        ret += sizeof(_clusterData.gridSize);
        ret += _clusterData.atomicVolume.required_allocation_size();
        ret += _clusterData.temperature.required_allocation_size(
            _clusterData.temperature.extent(0));
        ret += _clusterData.reactionRadius.required_allocation_size(
            _clusterData.reactionRadius.extent(0));
        ret += _clusterData.formationEnergy.required_allocation_size(
            _clusterData.formationEnergy.extent(0));
        ret += _clusterData.migrationEnergy.required_allocation_size(
            _clusterData.migrationEnergy.extent(0));
        ret += _clusterData.diffusionFactor.required_allocation_size(
            _clusterData.diffusionFactor.extent(0));
        ret += _clusterData.diffusionCoefficient.required_allocation_size(
            _clusterData.diffusionCoefficient.extent(0),
            _clusterData.diffusionCoefficient.extent(1));

        ret += sizeof(_reactionData.coeffExtent);
        ret += sizeof(_reactionData.numReactions);
        ret += _reactionData.productionCoeffs.required_allocation_size(
            _reactionData.productionCoeffs.extent(0),
            _reactionData.productionCoeffs.extent(1),
            _reactionData.productionCoeffs.extent(2),
            _reactionData.productionCoeffs.extent(3),
            _reactionData.productionCoeffs.extent(4));
        ret += _reactionData.dissociationCoeffs.required_allocation_size(
            _reactionData.dissociationCoeffs.extent(0),
            _reactionData.dissociationCoeffs.extent(1),
            _reactionData.dissociationCoeffs.extent(2),
            _reactionData.dissociationCoeffs.extent(3),
            _reactionData.dissociationCoeffs.extent(4));
        ret += _reactionData.rates.required_allocation_size(
            _reactionData.rates.extent(0), _reactionData.rates.extent(1));

        ret += _reactionData.connectivity.getDeviceMemorySize();

        ret += _reactions.required_allocation_size(_reactions.extent(0));

        return ret;
    }

    void
    syncClusterDataOnHost() override
    {
        _subpaving.syncTiles(plsm::onHost);
        auto mirror = ClusterDataMirror(_subpaving, this->_gridSize);
        Kokkos::deep_copy(mirror.atomicVolume, _clusterData.atomicVolume);
        Kokkos::deep_copy(mirror.temperature, _clusterData.temperature);
        Kokkos::deep_copy(mirror.momentIds, _clusterData.momentIds);
        Kokkos::deep_copy(mirror.reactionRadius, _clusterData.reactionRadius);
        Kokkos::deep_copy(mirror.formationEnergy, _clusterData.formationEnergy);
        Kokkos::deep_copy(mirror.migrationEnergy, _clusterData.migrationEnergy);
        Kokkos::deep_copy(mirror.diffusionFactor, _clusterData.diffusionFactor);
        Kokkos::deep_copy(mirror.diffusionCoefficient, _clusterData.diffusionCoefficient);
        _clusterDataMirror = mirror;
    }

    KOKKOS_INLINE_FUNCTION
    Cluster<plsm::OnDevice>
    findCluster(const Composition& comp, plsm::OnDevice context)
    {
        return Cluster<plsm::OnDevice>(_clusterData,
            _subpaving.findTileId(comp, context));
    }

    Cluster<plsm::OnHost>
    findCluster(const Composition& comp, plsm::OnHost context)
    {
        return Cluster<plsm::OnHost>(_clusterDataMirror,
            _subpaving.findTileId(comp, context));
    }

    KOKKOS_INLINE_FUNCTION
    auto
    findCluster(const Composition& comp)
    {
        return findCluster(comp, plsm::onDevice);
    }

    ClusterCommon<plsm::OnHost>
    getClusterCommon(std::size_t clusterId) const override
    {
        return ClusterCommon<plsm::OnHost>(_clusterDataMirror, clusterId);
    }

    ClusterCommon<plsm::OnHost>
    getSingleVacancy() override
    {
        Composition comp = Composition::zero();

        // Find the vacancy index
        constexpr auto speciesRangeNoI = getSpeciesRangeNoI();
        bool hasVacancy = false;
        Species vIndex;
        for (auto i : speciesRangeNoI) {
            if (isVacancy(i)) {
                hasVacancy = true;
                vIndex = i;
            }
        }

        // Update the composition if there is vacancy in the network
        if (hasVacancy) comp[vIndex] = 1;

        auto clusterId = findCluster(comp, plsm::onHost).getId();

        return ClusterCommon<plsm::OnHost>(_clusterDataMirror, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    Cluster<plsm::OnDevice>
    getCluster(std::size_t clusterId, plsm::OnDevice)
    {
        return Cluster<plsm::OnDevice>(_clusterData, clusterId);
    }

    Cluster<plsm::OnHost>
    getCluster(std::size_t clusterId, plsm::OnHost)
    {
        return Cluster<plsm::OnHost>(_clusterDataMirror, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getCluster(std::size_t clusterId)
    {
        return getCluster(clusterId, plsm::onDevice);
    }

    KOKKOS_INLINE_FUNCTION
    Subpaving&
    getSubpaving()
    {
        return _subpaving;
    }

    void
    computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex) override;

    void
    computeAllPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, std::size_t gridIndex) override;

    double
    getLargestRate() override;

    double
    getLeftSideRate(ConcentrationsView concentrations,
        std::size_t clusterId, std::size_t gridIndex) override;

    std::size_t
    getDiagonalFill(SparseFillMap& fillMap) override;

    /**
     * Get the total concentration of a given type of clusters.
     *
     * @param concentration The vector of concentrations
     * @param type The type of atom we want the concentration of
     * @param minSize The minimum number of atom to start counting
     * @return The total accumulated concentration times the size of the cluster
     */
    double
    getTotalConcentration(ConcentrationsView concentrations,
            Species type, AmountType minSize = 0);

    /**
     * Get the total concentration of a given type of clusters times their radius.
     *
     * @param concentration The vector of concentrations
     * @param type The type of atom we want the concentration of
     * @param minSize The minimum number of atom to start counting
     * @return The total accumulated concentration times the size of the cluster
     */
    double
    getTotalRadiusConcentration(ConcentrationsView concentrations,
            Species type, AmountType minSize = 0);

    /**
     * Get the total concentration of a given type of clusters times the number of atoms.
     *
     * @param concentration The vector of concentrations
     * @param type The type of atom we want the concentration of
     * @param minSize The minimum number of atom to start counting
     * @return The total accumulated concentration times the size of the cluster
     */
    double
    getTotalAtomConcentration(ConcentrationsView concentrations,
            Species type, AmountType minSize = 0);

    /**
     * Get the total concentration of a given type of clusters only if it is trapped in a vacancy.
     *
     * @param concentration The vector of concentrations
     * @param type The type of atom we want the concentration of
     * @param minSize The minimum number of atom to start counting
     * @return The total accumulated concentration times the size of the cluster
     */
    double
    getTotalTrappedAtomConcentration(ConcentrationsView concentrations,
            Species type, AmountType minSize = 0);

private:
    KOKKOS_INLINE_FUNCTION
    TImpl*
    asDerived()
    {
        return static_cast<TImpl*>(this);
    }

    void
    defineMomentIds();

    void
    generateClusterData(const ClusterGenerator& generator);

    void
    defineReactions();

    void
    updateDiffusionCoefficients();

    KOKKOS_INLINE_FUNCTION
    double
    getTemperature(std::size_t gridIndex) const noexcept
    {
        return _clusterData.temperature(gridIndex);
    }

private:
    Subpaving _subpaving;

    ClusterData _clusterData;
    ClusterDataMirror _clusterDataMirror;

    detail::ReactionData _reactionData;
    Kokkos::View<ReactionType*> _reactions;

    detail::ReactionNetworkWorker<TImpl> _worker;
};


namespace detail
{
template <typename TImpl>
struct ReactionNetworkWorker
{
    using Network = ReactionNetwork<TImpl>;
    using Types = ReactionNetworkTypes<TImpl>;
    using ClusterData = typename Types::ClusterData;
    using ClusterDataRef = typename Types::ClusterDataRef;
    using IndexType = typename Types::IndexType;

    Network& _nw;

    ReactionNetworkWorker(Network& network)
        :
        _nw(network)
    {
    }

    void
    updateDiffusionCoefficients();

    void
    generateClusterData(const typename Network::ClusterGenerator& generator);

    void
    defineMomentIds();

    struct ClusterSetsViewPair
    {
        using ClusterSet = typename Network::ReactionType::ClusterSet;
        Kokkos::View<ClusterSet*> prodClusterSets;
        Kokkos::View<ClusterSet*> dissClusterSets;
    };

    void
    defineReactions();

    std::size_t
    getDiagonalFill(typename Network::SparseFillMap& fillMap);
};
}
}
}

#include <experimental/ReactionGenerator.h>
#include <experimental/ReactionNetwork.inl>
