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
    _gridSize(gridSize),
    _temperature("Temperature", _gridSize),
    _numClusters(_subpaving.getNumberOfTiles(plsm::onDevice)),
    _diffusionCoefficient("Diffusion Coefficient", _numClusters, _gridSize),
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
            ////
            subpaving.syncAll(plsm::onHost);
            ////
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
ReactionNetwork<TImpl>::setTemperatures(const std::vector<double>& gridTemps)
{
    Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> tempsHost(
        gridTemps.data(), _gridSize);
    _temperature = Kokkos::View<double*>("Grid Temperatures", _gridSize);
    Kokkos::deep_copy(_temperature, tempsHost);

    updateDiffusionCoefficients();

    Kokkos::parallel_for(_reactions.extent(0), KOKKOS_LAMBDA (std::size_t i) {
        _reactions(i).updateRates();
    });
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::updateDiffusionCoefficients()
{
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Kokkos::parallel_for(Range2D({0, 0}, {_numClusters, _gridSize}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        _diffusionCoefficient(i,j) = _diffusionFactor(i) * std::exp(
            -_migrationEnergy(i) / (kBoltzmann * _temperature(j)));
    });
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
    
    // Set DOF
    _numDOFs = current - 1;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineReactions()
{
    _subpaving.syncTiles(plsm::onHost);
    auto tiles = _subpaving.getTiles(plsm::onHost);
    auto diffusionFactor = Kokkos::create_mirror_view(_diffusionFactor);
    Kokkos::deep_copy(diffusionFactor, _diffusionFactor);

    using ClusterSet = typename ReactionType::ClusterSet;
    ClusterSetsPair clusterSetsPair =
        asDerived()->defineReactionClusterSets(tiles, diffusionFactor);
    auto& prodClusterSets = clusterSetsPair.prodClusterSets;
    auto& dissClusterSets = clusterSetsPair.dissClusterSets;
    auto nProdReactions = prodClusterSets.size();
    auto nDissReactions = dissClusterSets.size();

    auto prodSetsHostView = 
        Kokkos::View<ClusterSet*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
            prodClusterSets.data(), nProdReactions);
    auto prodSets = Kokkos::View<ClusterSet*>("Production Cluster Sets",
        nProdReactions);
    Kokkos::deep_copy(prodSets, prodSetsHostView);
    auto dissSetsHostView =
        Kokkos::View<ClusterSet*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
            dissClusterSets.data(), nDissReactions);
    auto dissSets = Kokkos::View<ClusterSet*>("Dissociation Cluster Sets",
        nDissReactions);
    Kokkos::deep_copy(dissSets, dissSetsHostView);

    static constexpr auto cExt = getNumberOfSpeciesNoI() + 1;
    Kokkos::resize(_productionCoeffs, nProdReactions, cExt, cExt, 4, cExt);
    Kokkos::resize(_dissociationCoeffs, nDissReactions, cExt, 1, 3, cExt);

    auto numReactions = nProdReactions + nDissReactions;
    Kokkos::resize(_reactions, numReactions);
    _reactions = Kokkos::View<ReactionType*>("Reactions", numReactions);
    Kokkos::resize(_reactionRates, numReactions);

    using RType = typename ReactionType::Type;
    Kokkos::parallel_for(nProdReactions, KOKKOS_LAMBDA (std::size_t i) {
        const auto& clSet = prodSets(i);
        _reactions(i) = ReactionType(*this, i, RType::production,
            clSet.cluster0, clSet.cluster1, clSet.cluster2, clSet.cluster3);
    });
    Kokkos::parallel_for(nDissReactions, KOKKOS_LAMBDA (std::size_t i) {
        const auto& clSet = dissSets(i);
        auto id = i + nProdReactions;
        _reactions(id) = ReactionType(*this, id, RType::dissociation,
            clSet.cluster0, clSet.cluster1, clSet.cluster2, clSet.cluster3);
    });
}

template <typename TImpl>
std::size_t
ReactionNetwork<TImpl>::getDiagonalFill(SparseFillMap& fillMap)
{
    // Create the connectivity matrix initialized to invalid
    using Kokkos::ViewAllocateWithoutInitializing;
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Kokkos::View<std::size_t**> connectivity(
        ViewAllocateWithoutInitializing("connectivity"), _numDOFs, _numDOFs);
    Kokkos::parallel_for(Range2D({0, 0}, {_numDOFs, _numDOFs}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        connectivity(i,j) = invalid;
    });
    // Loop on each reaction to add its contribution to the connectivity matrix
    const auto& nReactions = _reactions.extent(0);
    Kokkos::parallel_for(nReactions, KOKKOS_LAMBDA (const std::size_t i) {
        _reactions(i).contributeConnectivity(connectivity);
    });

    // Initialize the inverse map to invalid as well
    _inverseMap = Kokkos::View<std::size_t**, Kokkos::MemoryUnmanaged>(
        ViewAllocateWithoutInitializing("_inverseMap"), _numDOFs, _numDOFs);
    Kokkos::parallel_for(Range2D({0, 0}, {_numDOFs, _numDOFs}),
            KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
        _inverseMap(i,j) = invalid;
    });
    // Transfer to fillMap, fill the inverse map,
    // and count the total number of partials
    std::size_t nPartials = 0;
    for (std::size_t i = 0; i < _numDOFs; ++i) {
        // Create a vector for ids
        std::vector<int> current;
        // Loop on this row
        for (std::size_t j = 0; j < _numDOFs; ++j)
        {
            if (connectivity(i,j) == invalid) {
                // This is the end of the row
                break;
            }
            // Add the value to the vector
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
