#pragma once

#include <boost/timer/timer.hpp>

#include <Constants.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(const Subpaving& subpaving,
        std::size_t gridSize, const IOptions& options)
    :
    _subpaving(subpaving),
    _gridSize(gridSize),
    _clusterData(_subpaving, _gridSize),
    _worker(*this)
{
    setInterstitialBias(options.getBiasFactor());
    setImpurityRadius(options.getImpurityRadius());
    setLatticeParameter(options.getLatticeParameter());

    {
    boost::timer::auto_cpu_timer t;
    generateClusterData(ClusterGenerator{options});
    defineMomentIds();
    }

    {
    boost::timer::auto_cpu_timer t;
    defineReactions();
    }
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(const Subpaving& subpaving,
        std::size_t gridSize)
    :
    ReactionNetwork(subpaving, gridSize, Options{})
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
        const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        std::size_t gridSize, const IOptions& options)
    :
    ReactionNetwork(
        [&]() -> Subpaving
        {
            boost::timer::auto_cpu_timer t;
            Region latticeRegion{};
            for (std::size_t i = 0; i < getNumberOfSpecies(); ++i) {
                latticeRegion[i] = Ival{0, maxSpeciesAmounts[i] + 1};
            }
            Subpaving sp(latticeRegion, subdivisionRatios);
            sp.refine(ClusterGenerator{options});
            return sp;
        }(),
        gridSize, options)
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
        const std::vector<AmountType>& maxSpeciesAmounts,
        std::size_t gridSize, const IOptions& options)
    :
    ReactionNetwork(maxSpeciesAmounts,
        [&maxSpeciesAmounts]() -> std::vector<SubdivisionRatio>
        {
            SubdivisionRatio ratio;
            for (std::size_t i = 0; i < getNumberOfSpecies(); ++i) {
                ratio[i] = maxSpeciesAmounts[i] + 1;
            }
            return {ratio};
        }(),
        gridSize, options)
{
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setLatticeParameter(double latticeParameter)
{
    _latticeParameter = asDerived()->checkLatticeParameter(latticeParameter);

    _atomicVolume =
        0.5 * _latticeParameter * _latticeParameter * _latticeParameter;
    auto mirror = Kokkos::create_mirror_view(_clusterData.atomicVolume);
    mirror(0) = _atomicVolume;
    Kokkos::deep_copy(_clusterData.atomicVolume, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setTemperatures(const std::vector<double>& gridTemps)
{
    Kokkos::View<const double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
        tempsHost(gridTemps.data(), _gridSize);
    Kokkos::deep_copy(_clusterData.temperature, tempsHost);

    updateDiffusionCoefficients();

    Kokkos::parallel_for(_reactions.extent(0), KOKKOS_LAMBDA (std::size_t i) {
        _reactions(i).updateRates();
    });
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::updateDiffusionCoefficients()
{
    _worker.updateDiffusionCoefficients();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::generateClusterData(const ClusterGenerator& generator)
{
    _worker.generateClusterData(generator);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::computeAllFluxes(ConcentrationsView concentrations,
    FluxesView fluxes, std::size_t gridIndex)
{
    // Get the extent of the reactions view
    const auto& nReactions = _reactions.extent(0);
    // Loop on the reactions
    Kokkos::parallel_for(nReactions, KOKKOS_LAMBDA (const std::size_t i) {
        _reactions(i).contributeFlux(concentrations, fluxes, gridIndex);
    });
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::computeAllPartials(ConcentrationsView concentrations,
    Kokkos::View<double*> values, std::size_t gridIndex)
{
    // Reset the values
    const auto& nValues = values.extent(0);
    // Loop on the reactions
    Kokkos::parallel_for(nValues, KOKKOS_LAMBDA (const std::size_t i) {
        values(i) = 0.0;
    });

    // Get the extent of the reactions view
    const auto& nReactions = _reactions.extent(0);
    // Loop on the reactions
    Kokkos::parallel_for(nReactions, KOKKOS_LAMBDA (const std::size_t i) {
        _reactions(i).contributePartialDerivatives(concentrations, values, gridIndex);
    });
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineMomentIds()
{
    _worker.defineMomentIds(_numDOFs);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineReactions()
{
    _worker.defineReactions();
}

template <typename TImpl>
std::size_t
ReactionNetwork<TImpl>::getDiagonalFill(SparseFillMap& fillMap)
{
    return _worker.getDiagonalFill(fillMap);
}

namespace detail
{
template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::updateDiffusionCoefficients()
{
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    auto data = _nw._clusterData;
    Kokkos::parallel_for(Range2D({0, 0}, {data.numClusters, data.gridSize}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        data.diffusionCoefficient(i,j) = data.diffusionFactor(i) * std::exp(
            -data.migrationEnergy(i) / (kBoltzmann * data.temperature(j)));
    });
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::generateClusterData(
    const typename Network::ClusterGenerator& generator)
{
    auto nClusters = _nw._clusterData.numClusters;
    _nw._clusterData.formationEnergy = Kokkos::View<double*>("Formation Energy", nClusters);
    _nw._clusterData.migrationEnergy = Kokkos::View<double*>("Migration Energy", nClusters);
    _nw._clusterData.diffusionFactor = Kokkos::View<double*>("Diffusion Factor", nClusters);
    _nw._clusterData.reactionRadius = Kokkos::View<double*>("Reaction Radius", nClusters);

    ClusterData data(_nw._clusterData);
    auto latticeParameter = _nw.getLatticeParameter();
    auto interstitialBias = _nw.getInterstitialBias();
    auto impurityRadius = _nw.getImpurityRadius();
    Kokkos::parallel_for(nClusters, KOKKOS_LAMBDA (const std::size_t i) {
        auto cluster = data.getCluster(i);
        data.formationEnergy(i) = generator.getFormationEnergy(cluster);
        data.migrationEnergy(i) = generator.getMigrationEnergy(cluster);
        data.diffusionFactor(i) = generator.getDiffusionFactor(cluster);
        data.reactionRadius(i) = generator.getReactionRadius(cluster,
            latticeParameter, interstitialBias, impurityRadius);
    });
    Kokkos::fence();
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::defineMomentIds(std::size_t& numDOFs)
{
    constexpr auto speciesRangeNoI = Network::getSpeciesRangeNoI();

    ClusterDataRef data(_nw._clusterData);

    auto nClusters = data.numClusters;
    auto counts = Kokkos::View<std::size_t*>("Moment Id Counts", nClusters);

    std::size_t nMomentIds = 0;
    Kokkos::parallel_reduce(nClusters,
        KOKKOS_LAMBDA (const std::size_t i, std::size_t& running)
        {
            const auto& reg = data.getCluster(i).getRegion();
            std::size_t count = 0;
            for (auto k : speciesRangeNoI) {
                if (reg[k].length() != 1) {
                    ++count;
                }
            }
            running += count;
            counts(i) = count;
        },
        nMomentIds);

    Kokkos::parallel_scan(nClusters, ExclusiveScanFunctor{counts});

    Kokkos::parallel_for(nClusters, KOKKOS_LAMBDA (const std::size_t i) {
        const auto& reg = data.getCluster(i).getRegion();
        std::size_t current = counts(i);
        for (auto k : speciesRangeNoI) {
            if (reg[k].length() == 1) {
                data.momentIds(i, k()) = Network::invalid;
            }
            else {
                data.momentIds(i, k()) = current;
                ++current;
            }
        }
    });

    Kokkos::fence();
    numDOFs = nClusters + nMomentIds;
}

template <typename TImpl>
typename ReactionNetworkWorker<TImpl>::ClusterSetsViewPair
ReactionNetworkWorker<TImpl>::defineReactionClusterSets()
{
    using ClusterSet = typename ClusterSetsViewPair::ClusterSet;
    ClusterSetsViewPair ret;

    auto subpaving = _nw._subpaving;
    auto tiles = subpaving.getTiles(plsm::onDevice);
    std::size_t numClusters = tiles.extent(0);

    detail::UpperTriangle<ClusterSet> prodSet("Temp Production Set",
        numClusters);
    auto cap = prodSet.size();
    detail::UpperTriangle<ClusterSet> dissSet("Temp Dissociation Set",
        numClusters);

    auto diffusionFactor = _nw._clusterData.diffusionFactor;
    auto validator = _nw.asDerived()->getReactionValidator();
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Kokkos::parallel_for(Range2D({0, 0}, {numClusters, numClusters}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        if (j < i) {
            return;
        }
        if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
            return;
        }
        validator(i, j, subpaving, prodSet, dissSet);
    });
    Kokkos::fence();

    auto ids = Kokkos::View<std::size_t*>("Reaction Ids", cap);
    std::size_t numProdReactions = 0;
    Kokkos::parallel_scan(cap,
            KOKKOS_LAMBDA (std::size_t i, std::size_t& update,
                const bool finalPass) {
        if (prodSet(i).valid()) {
            if (finalPass) {
                ids(i) = update;
            }
            update += 1;
        }
    }, numProdReactions);
    auto prodSetsView = Kokkos::View<ClusterSet*>("Production Cluster Sets",
        numProdReactions);
    Kokkos::parallel_for(cap, KOKKOS_LAMBDA (std::size_t i) {
        if (prodSet(i).valid()) {
            prodSetsView(ids(i)) = prodSet(i);
        }
    });

    std::size_t numDissReactions = 0;
    Kokkos::parallel_scan(cap,
            KOKKOS_LAMBDA (std::size_t i, std::size_t& update,
                const bool finalPass) {
        if (dissSet(i).valid()) {
            if (finalPass) {
                ids(i) = update;
            }
            update += 1;
        }
    }, numDissReactions);
    auto dissSetsView = Kokkos::View<ClusterSet*>("Dissociation Cluster Sets",
        numDissReactions);
    Kokkos::parallel_for(cap, KOKKOS_LAMBDA (std::size_t i) {
        if (dissSet(i).valid()) {
            dissSetsView(ids(i)) = dissSet(i);
        }
    });

    return {prodSetsView, dissSetsView};
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::defineReactions()
{
    using ReactionType = typename Network::ReactionType;
    using ClusterSet = typename ReactionType::ClusterSet;

    auto clusterSetsPair = defineReactionClusterSets();
    auto prodSets = clusterSetsPair.prodClusterSets;
    auto nProdReactions = prodSets.extent(0);
    auto dissSets = clusterSetsPair.dissClusterSets;
    auto nDissReactions = dissSets.extent(0);

    static constexpr auto cExt = Network::getNumberOfSpeciesNoI() + 1;
    _nw._productionCoeffs = Kokkos::View<double*****>("Production Coefficients",
        nProdReactions, cExt, cExt, 4, cExt);
    _nw._dissociationCoeffs = Kokkos::View<double*****>(
        "Dissociation Coefficients", nDissReactions, cExt, 1, 3, cExt);

    auto numReactions = nProdReactions + nDissReactions;
    _nw._reactionRates = Kokkos::View<double**>("Reaction Rates", numReactions,
        _nw._gridSize);
    _nw._reactions = Kokkos::View<ReactionType*>("Reactions", numReactions);

    // Initialize the inverse map to invalid
    _nw._inverseMap = Kokkos::View<std::size_t**>(
        Kokkos::ViewAllocateWithoutInitializing("_inverseMap"),
        _nw._numDOFs, _nw._numDOFs);
    auto invMap = _nw._inverseMap;
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Kokkos::parallel_for(Range2D({0, 0}, {_nw._numDOFs, _nw._numDOFs}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        invMap(i,j) = Network::invalid;
    });

    auto reactions = _nw._reactions;
    ReactionDataRef reactionData{_nw._productionCoeffs, _nw._dissociationCoeffs,
        _nw._reactionRates, _nw._inverseMap};

    using RType = typename ReactionType::Type;
    auto data = _nw._clusterData;
    Kokkos::parallel_for(nProdReactions, KOKKOS_LAMBDA (std::size_t i) {
        const auto& clSet = prodSets(i);
        reactions(i) = ReactionType(reactionData, data, i, RType::production,
            clSet.cluster0, clSet.cluster1, clSet.cluster2, clSet.cluster3);
    });
    Kokkos::parallel_for(nDissReactions, KOKKOS_LAMBDA (std::size_t i) {
        const auto& clSet = dissSets(i);
        auto id = i + nProdReactions;
        reactions(id) = ReactionType(reactionData, data, id, RType::dissociation,
            clSet.cluster0, clSet.cluster1, clSet.cluster2, clSet.cluster3);
    });
    Kokkos::fence();
}

template <typename TImpl>
std::size_t
ReactionNetworkWorker<TImpl>::getDiagonalFill(
    typename Network::SparseFillMap& fillMap)
{
    // Create the connectivity matrix initialized to invalid
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Kokkos::View<std::size_t**> connectivity(
        Kokkos::ViewAllocateWithoutInitializing("connectivity"),
        _nw._numDOFs, _nw._numDOFs);
    Kokkos::parallel_for(Range2D({0, 0}, {_nw._numDOFs, _nw._numDOFs}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        connectivity(i,j) = Network::invalid;
    });
    // Loop on each reaction to add its contribution to the connectivity matrix
    const auto& nReactions = _nw._reactions.extent(0);
    Kokkos::parallel_for(nReactions, KOKKOS_LAMBDA (const std::size_t i) {
        _nw._reactions(i).contributeConnectivity(connectivity);
    });

    // Transfer to fillMap, fill the inverse map,
    // and count the total number of partials
    std::size_t nPartials = 0;
    for (std::size_t i = 0; i < _nw._numDOFs; ++i) {
        // Create a vector for ids
        std::vector<int> current;
        // Loop on this row
        for (std::size_t j = 0; j < _nw._numDOFs; ++j)
        {
            if (connectivity(i,j) == Network::invalid) {
                // This is the end of the row
                break;
            }
            // Add the value to the vector
            current.push_back((int) connectivity(i,j));
            // Update the inverse map
            _nw._inverseMap(i,connectivity(i,j)) = nPartials;
            
            // Count
            nPartials++;
        }
        // Add the current vector to fillMap
        fillMap[(int) i] = current;
    }

    return nPartials;
}
}
}
}
