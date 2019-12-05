#pragma once

#include <Constants.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(Subpaving&& subpaving,
        std::size_t gridSize, const IOptions& options)
    :
    _subpaving(std::move(subpaving)),
    _temperature("Temperature", gridSize),
    _numClusters(_subpaving.getNumberOfTiles(plsm::onDevice)),
    _diffusionCoefficient("Diffusion Coefficient", _numClusters, gridSize),
    _reactionRadius("Reaction Radius", _numClusters),
    _formationEnergy("Formation Energy", _numClusters),
    _migrationEnergy("Migration Energy", _numClusters),
    _diffusionFactor("Diffusion Factor", _numClusters)
{
    setInterstitialBias(options.getBiasFactor());
    setImpurityRadius(options.getImpurityRadius());
    setLatticeParameter(options.getLatticeParameter());

    generateClusterData(ClusterGenerator{options});
    defineMomentIds();

    defineReactions();
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(Subpaving&& subpaving,
        std::size_t gridSize)
    :
    ReactionNetwork(std::move(subpaving), gridSize, Options{})
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
        const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        std::size_t gridSize, const IOptions& options)
    :
    ReactionNetwork(
        [&maxSpeciesAmounts, &subdivisionRatios, &options]()
        {
            Region latticeRegion;
            for (std::size_t i = 0; i < getNumberOfSpecies(); ++i) {
                latticeRegion[i] = Ival{0, maxSpeciesAmounts[i] + 1};
            }
            Subpaving subpaving(latticeRegion, subdivisionRatios);
            subpaving.refine(ClusterGenerator{options});
            return subpaving;
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
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::generateClusterData(const ClusterGenerator& generator)
{
    Kokkos::parallel_for(_numClusters, KOKKOS_LAMBDA (const std::size_t i) {
        auto cluster = getCluster(i);
        _formationEnergy(i) = generator.getFormationEnergy(cluster);
        _migrationEnergy(i) = generator.getMigrationEnergy(cluster);
        _diffusionFactor(i) = generator.getDiffusionFactor(cluster);
    });
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
    constexpr auto invalid = plsm::invalid<std::size_t>;

    _subpaving.syncAll(plsm::onHost);
    auto tiles = _subpaving.getTiles(plsm::onHost);
    auto nClusters = _subpaving.getNumberOfTiles(plsm::onHost);

    //FIXME: _momentIds data lives on the device this way
    _momentIds = Kokkos::View<std::size_t*[4]>("Moment Ids", nClusters);

    auto current = nClusters;
    for (std::size_t c = 0; c < nClusters; ++c) {
        const auto& reg = tiles(c).getRegion();
        for (auto k : getSpeciesRangeNoI()) {
            if (reg[k].length() == 1) {
                _momentIds(c, k()) = invalid;
            }
            else {
                _momentIds(c, k()) = current;
                ++current;
            }
        }
    }
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineReactions()
{
    _subpaving.syncTiles(plsm::onHost);
    auto tiles = _subpaving.getTiles(plsm::onHost);
    auto diffusionFactor = Kokkos::create_mirror_view(_diffusionFactor);
    Kokkos::deep_copy(diffusionFactor, _diffusionFactor);

    using ClusterAssoc = typename ReactionType::ClusterAssoc;
    std::vector<ClusterAssoc> clusterSets =
        asDerived()->defineReactionClusterSets(tiles, diffusionFactor);

    auto numReactions = clusterSets.size();
    typename Kokkos::View<ClusterAssoc*, Kokkos::MemoryUnmanaged>::HostMirror
        clSetsHostView(clusterSets.data(), numReactions);
    Kokkos::View<ClusterAssoc*> clSets("Cluster Sets", numReactions);
    Kokkos::deep_copy(clSets, clSetsHostView);
    _reactions = Kokkos::View<ReactionType*>("Reactions", numReactions);

    Kokkos::parallel_for(numReactions, KOKKOS_LAMBDA (std::size_t i) {
        const auto& clSet = clSets(i);
        _reactions(i) = ReactionType(*this, i, clSet.reactionType,
            clSet.cluster0, clSet.cluster1, clSet.cluster2, clSet.cluster3);
    });
}

template <typename TImpl>
size_t
ReactionNetwork<TImpl>::getDiagonalFill(SparseFillMap& fillMap)
{
    // TODO: initialize connectivity to invalid
    auto connectivity = Kokkos::View<std::size_t**, Kokkos::MemoryUnmanaged>(_numDOFs, _numDOFs);
    // Loop on each reaction to add its contribution to the connectivity matrix
    const auto& nReactions = _reactions.extent(0);
    Kokkos::parallel_for(nReactions, KOKKOS_LAMBDA (const std::size_t i) {
        _reactions(i).contributeConnectivity(connectivity);
    });

    // Transfer to fillMap, fill the inverse map,
    // and count the total number of partials
    // TODO: should it be initialized to invalid as well?
    _inverseMap = Kokkos::View<std::size_t**, Kokkos::MemoryUnmanaged>(_numDOFs, _numDOFs);
    std::size_t nPartials = 0;
    for (std::size_t i = 0; i < _numDOFs; ++i) 
    {
        // Create a vector for ids
        std::vector<int> current;
        // Loop on this row
        for (std::size_t j = 0; j < _numDOFs; ++j)
        {
            if (connectivity(i,j) == invalid) 
            {
                // This is the end of the row
                break;
            }
            // Add the value to the vector
            // TODO: Should we use emplace_back instead?
            current.push_back((int) connectivity(i,j));
            // Update the inverse map
            _inverseMap(i,connectivity(i,j)) = nPartials;
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
