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
    _diffusionFactor("Diffusion Factor", _numClusters),
    _reactions("Reactions", 1000)
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

// TODO
// template <typename TImpl>
// ReactionNetwork<TImpl>::ReactionNetwork(< subpaving region and subdivision infos >)
// {
// }

template <typename TImpl>
void
ReactionNetwork<TImpl>::setLatticeParameter(double latticeParameter)
{
    _latticeParameter = latticeParameter;
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
    //TODO
}

template <typename TImpl>
size_t
ReactionNetwork<TImpl>::getDiagonalFill(SparseFillMap& fillMap)
{
    //TODO
    return 0;
}
}
}
