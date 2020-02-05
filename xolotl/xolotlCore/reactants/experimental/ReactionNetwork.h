#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

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

template <typename TData>
class UpperTriangle
{
public:
    UpperTriangle(const std::string& label, std::size_t N)
        :
        _N(N),
        _data(label, ((N+1)*N)/2)
    {
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    size() const noexcept
    {
        return _data.extent(0);
    }

    KOKKOS_INLINE_FUNCTION
    TData&
    operator()(std::size_t i, std::size_t j) const
    {
        assert(j >= i);
        assert(i < _N);
        assert(j < _N);
        auto id = i*_N + j - ((i+1)*i)/2;
        return _data(id);
    }

    KOKKOS_INLINE_FUNCTION
    TData&
    operator()(std::size_t i) const
    {
        assert(i < size());
        return _data(i);
    }

private:
    std::size_t _N;
    Kokkos::View<TData*> _data;
};
}

template <typename TImpl>
class ReactionNetwork : public IReactionNetwork
{
    friend class detail::ReactionNetworkWorker<TImpl>;

    static constexpr auto invalid = plsm::invalid<std::size_t>;

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
    using ConnectivityView = typename IReactionNetwork::ConnectivityView;
    using SparseFillMap = typename IReactionNetwork::SparseFillMap;
    using ClusterData = typename Types::ClusterData;
    using ClusterDataMirror = typename Types::ClusterDataMirror;
    using ClusterDataRef = typename Types::ClusterDataRef;
    using InverseMap = Kokkos::View<std::size_t**>;

    template <typename PlsmContext>
    using Cluster = Cluster<Subpaving, PlsmContext>;

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
    setTemperatures(const std::vector<double>& gridTemperatures) override;

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

    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getReactionRates(std::size_t reactionId)
    {
        return Kokkos::subview(_reactionRates, reactionId, Kokkos::ALL);
    }

    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getReactionCoefficients(std::size_t reactionId)
    {
        if (reactionId < _productionCoeffs.extent(0)) {
            return Kokkos::subview(_productionCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
        else {
            reactionId -= _productionCoeffs.extent(0);
            return Kokkos::subview(_dissociationCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
    }

    void
    computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex) override;

    void
    computeAllPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, std::size_t gridIndex) override;

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
    getTotalAtomConcentration(ConcentrationsView concentrations,
    		Species type, std::size_t minSize = 0);

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
    		Species type, std::size_t minSize = 0);

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

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getReactionIndex(std::size_t rowId, std::size_t colId) const noexcept
    {
        return _inverseMap(rowId, colId);
    }

private:
    Subpaving _subpaving;

    ClusterData _clusterData;
    ClusterDataMirror _clusterDataMirror;

    Kokkos::View<ReactionType*> _reactions;
    Kokkos::View<double*****> _productionCoeffs;
    Kokkos::View<double*****> _dissociationCoeffs;
    Kokkos::View<double**> _reactionRates;

    // TODO: the original code uses an actual map here because it is sparse
    //       Look into Kokkos::Crs
    Kokkos::View<std::size_t**> _inverseMap;

    detail::ReactionNetworkWorker<TImpl> _worker;
};


namespace detail
{
template <typename TImpl>
struct ReactionNetworkWorker
{
    class ExclusiveScanFunctor
    {
    public:
        ExclusiveScanFunctor(Kokkos::View<std::size_t*> data)
            :
            _data(data)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(std::size_t index, std::size_t& update, const bool finalPass)
            const
        {
            const auto temp = _data(index);
            if (finalPass) {
                _data(index) = update;
            }
            update += temp;
        }

        KOKKOS_INLINE_FUNCTION
        void
        init(std::size_t& update) const
        {
          update = 0;
        }

        KOKKOS_INLINE_FUNCTION
        void
        join(volatile std::size_t& update, std::size_t input) const
        {
          update += input;
        }

    private:
        Kokkos::View<std::size_t*> _data;
    };

    using Network = ReactionNetwork<TImpl>;
    using ClusterData = typename Network::ClusterData;
    using ClusterDataRef = typename Network::ClusterDataRef;

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
    defineMomentIds(std::size_t& numDOFs);

    struct ClusterSetsViewPair
    {
        using ClusterSet = typename Network::ReactionType::ClusterSet;
        Kokkos::View<ClusterSet*> prodClusterSets;
        Kokkos::View<ClusterSet*> dissClusterSets;
    };

    ClusterSetsViewPair
    defineReactionClusterSets();

    void
    defineReactions();

    std::size_t
    getDiagonalFill(typename Network::SparseFillMap& fillMap);
};
}
}
}

#include <experimental/ReactionNetwork.inl>
