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
#include <experimental/ReactionCollection.h>
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
    using ClusterGenerator = typename Traits::ClusterGenerator;
    using AmountType = typename IReactionNetwork::AmountType;
    using IndexType = typename IReactionNetwork::IndexType;
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
    using ProductionReactionType = typename Traits::ProductionReactionType;
    using DissociationReactionType = typename Traits::DissociationReactionType;
    using SinkReactionType = typename Traits::SinkReactionType;
    using ReSolutionReactionType = typename Traits::ReSolutionReactionType;
    using ReactionCollection = detail::ReactionCollection<Traits>;

    template <typename PlsmContext>
    using Cluster = Cluster<TImpl, PlsmContext>;

    ReactionNetwork() = default;

    ReactionNetwork(const Subpaving& subpaving, IndexType gridSize,
        const IOptions& options);

    ReactionNetwork(const Subpaving& subpaving, IndexType gridSize);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        IndexType gridSize, const IOptions& options);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        IndexType gridSize, const IOptions& options);

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
    setFissionRate(double rate) override;

    void
    setIsReaction(bool reaction) override;

    void
    setIsReSolution(bool reso) override;

    void
    setGridSize(IndexType gridSize) override;

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
        ret += _reactions.getDeviceMemorySize();

        return ret;
    }

    void
    syncClusterDataOnHost() override
    {
        _subpaving.syncTiles(plsm::onHost);
        auto mirror = ClusterDataMirror(_subpaving, this->_gridSize);
        Kokkos::deep_copy(mirror.atomicVolume, _clusterData.atomicVolume);
        Kokkos::deep_copy(mirror.latticeParameter, _clusterData.latticeParameter);
        Kokkos::deep_copy(mirror.fissionRate, _clusterData.fissionRate);
        Kokkos::deep_copy(mirror.isReaction, _clusterData.isReaction);
        Kokkos::deep_copy(mirror.isReSolution, _clusterData.isReSolution);
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
    getClusterCommon(IndexType clusterId) const override
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

    Bounds
    getAllClusterBounds() override
    {
        // Create the object to return
        Bounds bounds;

        // Loop on all the clusters
        constexpr auto speciesRange = getSpeciesRange();
        auto tiles = _subpaving.getTiles(plsm::onHost);
        for (IndexType i = 0; i < _numClusters; ++i) {
            const auto& clReg = tiles(i).getRegion();
            Composition lo = clReg.getOrigin();
            Composition hi = clReg.getUpperLimitPoint();
            std::vector<AmountType> boundVector;
            for (auto j : speciesRange) {
                boundVector.push_back(lo[j]);
                boundVector.push_back(hi[j] - 1);
            }
            bounds.push_back(boundVector);
        }
        return bounds;
    }

    PhaseSpace
    getPhaseSpace() override
    {
        // Create the object to return
        PhaseSpace space;

        // Loop on all the clusters
        constexpr auto speciesRange = getSpeciesRange();
        for (auto j : speciesRange) {
            space.push_back(toString(j));
        }
        return space;
    }


    KOKKOS_INLINE_FUNCTION
    Cluster<plsm::OnDevice>
    getCluster(IndexType clusterId, plsm::OnDevice)
    {
        return Cluster<plsm::OnDevice>(_clusterData, clusterId);
    }

    Cluster<plsm::OnHost>
    getCluster(IndexType clusterId, plsm::OnHost)
    {
        return Cluster<plsm::OnHost>(_clusterDataMirror, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getCluster(IndexType clusterId)
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
        IndexType gridIndex) override;

    void
    computeAllPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, IndexType gridIndex) override;

    double
    getLargestRate() override;

    double
    getLeftSideRate(ConcentrationsView concentrations,
        IndexType clusterId, IndexType gridIndex) override;

    IndexType
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
    getTemperature(IndexType gridIndex) const noexcept
    {
        return _clusterData.temperature(gridIndex);
    }

private:
    Subpaving _subpaving;
    ClusterDataMirror _clusterDataMirror;

    ReactionCollection _reactions;
    Kokkos::View<SinkReactionType*> _sinkReactions;
    Kokkos::View<ReSolutionReactionType*> _resoReactions;

    detail::ReactionNetworkWorker<TImpl> _worker;

protected:
    detail::ReactionData _reactionData;

    ClusterData _clusterData;
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
    using ReactionCollection =
        detail::ReactionCollection<typename Types::Traits>;

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

    void
    defineReactions();

    IndexType
    getDiagonalFill(typename Network::SparseFillMap& fillMap);
};
}
}
}

#include <experimental/ReactionGenerator.h>
#include <experimental/ReactionNetwork.inl>
