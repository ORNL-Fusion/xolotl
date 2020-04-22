#pragma once

#include <boost/timer/timer.hpp>

#include <Constants.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(const Subpaving& subpaving,
        IndexType gridSize, const IOptions& options)
    :
    IReactionNetwork(gridSize),
    _subpaving(subpaving),
    _clusterData(_subpaving, gridSize),
    _worker(*this)
{
    // Set constants
    setInterstitialBias(options.getBiasFactor());
    setImpurityRadius(options.getImpurityRadius());
    setLatticeParameter(options.getLatticeParameter());
    setFissionRate(options.getFluxAmplitude());
    auto map = options.getProcesses();
    setEnableStdReaction(map["reaction"]);
    setEnableReSolution(map["resolution"]);
    
    auto tiles = subpaving.getTiles(plsm::onDevice);
    _numClusters = tiles.extent(0);
    
//    // PRINT ALL THE CLUSTERS
//    constexpr auto speciesRange = getSpeciesRange();
//    for (IndexType i = 0; i < _numClusters; ++i) {
//        const auto& clReg = tiles(i).getRegion();
//        Composition lo = clReg.getOrigin();
//        Composition hi = clReg.getUpperLimitPoint();
//
//        std::cout << i << ": " << std::endl;
//        for (auto j : speciesRange) std::cout << lo[j] << " ";
//        std::cout << std::endl;
//        for (auto j : speciesRange) std::cout << hi[j] - 1 << " ";
//        std::cout << std::endl;
//    }
    std::cout << "num: " << _numClusters << std::endl;

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
        IndexType gridSize)
    :
    ReactionNetwork(subpaving, gridSize, Options{})
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
        const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        IndexType gridSize, const IOptions& options)
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
        IndexType gridSize, const IOptions& options)
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
    auto lParam = asDerived()->checkLatticeParameter(latticeParameter);
    this->_latticeParameter = lParam;
    auto mirror = Kokkos::create_mirror_view(_clusterData.latticeParameter);
    mirror(0) = this->_latticeParameter;
    Kokkos::deep_copy(_clusterData.latticeParameter, mirror);

    this->_atomicVolume =
        0.5 * lParam * lParam * lParam;
    mirror = Kokkos::create_mirror_view(_clusterData.atomicVolume);
    mirror(0) = this->_atomicVolume;
    Kokkos::deep_copy(_clusterData.atomicVolume, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setFissionRate(double rate)
{
    this->_fissionRate = rate;
    auto mirror = Kokkos::create_mirror_view(_clusterData.fissionRate);
    mirror(0) = this->_fissionRate;
    Kokkos::deep_copy(_clusterData.fissionRate, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableStdReaction(bool reaction)
{
    this->_enableStdReaction = reaction;
    auto mirror = Kokkos::create_mirror_view(_clusterData.enableStdReaction);
    mirror(0) = this->_enableStdReaction;
    Kokkos::deep_copy(_clusterData.enableStdReaction, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableReSolution(bool reso)
{
    this->_enableReSolution = reso;
    auto mirror = Kokkos::create_mirror_view(_clusterData.enableReSolution);
    mirror(0) = this->_enableReSolution;
    Kokkos::deep_copy(_clusterData.enableReSolution, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setGridSize(IndexType gridSize)
{
    _gridSize = gridSize;
    _clusterData.setGridSize(_gridSize);
    _clusterDataMirror.setGridSize(_gridSize);
    _reactionData.setGridSize(_gridSize);

    auto reactionData = _reactionData;
    auto clusterData = _clusterData;
    _reactions.apply(DEVICE_LAMBDA (auto&& reaction) {
        reaction.updateData(reactionData, clusterData);
    });
    auto sinks = _sinkReactions;
    auto numSinks = sinks.extent(0);
    Kokkos::parallel_for(numSinks, KOKKOS_LAMBDA (const IndexType i) {
        sinks(i).updateData(reactionData, clusterData);
    });
    auto resos = _resoReactions;
    auto numReSos = resos.extent(0);
    Kokkos::parallel_for(numReSos, KOKKOS_LAMBDA (const IndexType i) {
        resos(i).updateData(reactionData, clusterData);
    });
    Kokkos::fence();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setTemperatures(const std::vector<double>& gridTemps)
{
    Kokkos::View<const double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
        tempsHost(gridTemps.data(), this->_gridSize);
    Kokkos::deep_copy(_clusterData.temperature, tempsHost);

    updateDiffusionCoefficients();

    _reactions.apply(DEVICE_LAMBDA (auto&& reaction) {
        reaction.updateRates();
    });
    auto sinks = _sinkReactions;
    auto numSinks = sinks.extent(0);
    Kokkos::parallel_for(numSinks, KOKKOS_LAMBDA (const IndexType i) {
        sinks(i).updateRates();
    });
    auto resos = _resoReactions;
    auto numReSos = resos.extent(0);
    Kokkos::parallel_for(numReSos, KOKKOS_LAMBDA (const IndexType i) {
        resos(i).updateRates();
    });
    Kokkos::fence();
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
    FluxesView fluxes, IndexType gridIndex)
{
    _reactions.apply(DEVICE_LAMBDA (auto&& reaction) {
        reaction.contributeFlux(concentrations, fluxes, gridIndex);
    });
    auto sinks = _sinkReactions;
    auto numSinks = sinks.extent(0);
    Kokkos::parallel_for(numSinks, KOKKOS_LAMBDA (const IndexType i) {
        sinks(i).contributeFlux(concentrations, fluxes, gridIndex);
    });
    auto resos = _resoReactions;
    auto numReSos = resos.extent(0);
    Kokkos::parallel_for(numReSos, KOKKOS_LAMBDA (const IndexType i) {
        resos(i).contributeFlux(concentrations, fluxes, gridIndex);
    });

    Kokkos::fence();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::computeAllPartials(ConcentrationsView concentrations,
    Kokkos::View<double*> values, IndexType gridIndex)
{
    // Reset the values
    const auto& nValues = values.extent(0);
    // Loop on the reactions
    Kokkos::parallel_for(nValues, KOKKOS_LAMBDA (const IndexType i) {
        values(i) = 0.0;
    });

    auto connectivity = _reactionData.connectivity;
    _reactions.apply(DEVICE_LAMBDA (auto&& reaction) {
        reaction.contributePartialDerivatives(concentrations, values,
            connectivity, gridIndex);
    });
    auto sinks = _sinkReactions;
    auto numSinks = sinks.extent(0);
    Kokkos::parallel_for(numSinks, KOKKOS_LAMBDA (const IndexType i) {
        sinks(i).contributePartialDerivatives(concentrations, values,
                connectivity, gridIndex);
    });
    auto resos = _resoReactions;
    auto numReSos = resos.extent(0);
    Kokkos::parallel_for(numReSos, KOKKOS_LAMBDA (const IndexType i) {
        resos(i).contributePartialDerivatives(concentrations, values,
                connectivity, gridIndex);
    });
    
    Kokkos::fence();
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getLargestRate()
{
    // Get the extent of the reactions rates view and grid
    auto reactionRates = _reactionData.rates;
    const auto& nRates = reactionRates.extent(0);
    double largestRate = 0.0;
    // Loop on all the rates to get the maximum
    using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    Kokkos::parallel_reduce(Range2D({0, 0}, {nRates, _gridSize}),
            KOKKOS_LAMBDA (IndexType i, IndexType j, double &max) {
        if (reactionRates(i, j) > max) max = reactionRates(i, j);
    }, Kokkos::Max<double>(largestRate));

    return largestRate;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getLeftSideRate(ConcentrationsView concentrations,
    IndexType clusterId, IndexType gridIndex)
{
    // Get the extent of the reactions
    double leftSideRate = 0.0;
    // Loop on all the rates to get the maximum
    _reactions.reduce(DEVICE_LAMBDA (auto&& reaction, double& lsum) {
        lsum += reaction.contributeLeftSideRate(concentrations, clusterId,
            gridIndex);
    }, leftSideRate);

    return leftSideRate;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalConcentration(ConcentrationsView concentrations,
        Species type, AmountType minSize)
{
    auto tiles = _subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    Kokkos::parallel_reduce(_numClusters,
            KOKKOS_LAMBDA (IndexType i, double &lsum) {
    	const Region& clReg = tiles(i).getRegion();
    	const auto factor = clReg.volume() / clReg[type].length();
    	for (AmountType j : makeIntervalRange(clReg[type])) {
    		if (j >= minSize) lsum += concentrations(i) * factor;
    	}
    }, conc);

    return conc;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalRadiusConcentration(ConcentrationsView concentrations,
        Species type, AmountType minSize)
{
    auto tiles = _subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    auto clusterData = _clusterData;
    Kokkos::parallel_reduce(_numClusters,
            KOKKOS_LAMBDA (IndexType i, double &lsum) {
    	const Region& clReg = tiles(i).getRegion();
    	const auto factor = clReg.volume() / clReg[type].length();
    	for (AmountType j : makeIntervalRange(clReg[type])) {
    		if (j >= minSize) lsum += concentrations(i) * clusterData.reactionRadius(i) * factor;
    	}
    }, conc);

    return conc;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalAtomConcentration(ConcentrationsView concentrations,
        Species type, AmountType minSize)
{
    auto tiles = _subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    Kokkos::parallel_reduce(_numClusters,
            KOKKOS_LAMBDA (IndexType i, double &lsum) {
    	const Region& clReg = tiles(i).getRegion();
    	const auto factor = clReg.volume() / clReg[type].length();
    	for (AmountType j : makeIntervalRange(clReg[type])) {
    		if (j >= minSize) lsum += concentrations(i) * j * factor;
    	}
    }, conc);

    return conc;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalTrappedAtomConcentration(ConcentrationsView concentrations,
        Species type, AmountType minSize)
{
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

    // Return 0 if there is not vacancy in the network
    if (!hasVacancy) return 0.0;

    auto tiles = _subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    Kokkos::parallel_reduce(_numClusters,
            KOKKOS_LAMBDA (IndexType i, double &lsum) {
        const Region& clReg = tiles(i).getRegion();
        if (clReg[vIndex].begin() > 0) {
            const auto factor = clReg.volume() / clReg[type].length();
            for (AmountType j : makeIntervalRange(clReg[type])) {
                if (j >= minSize) lsum += concentrations(i) * j * factor;
            }
        }
    }, conc);

    return conc;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineMomentIds()
{
    _worker.defineMomentIds();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineReactions()
{
    _worker.defineReactions();
}

template <typename TImpl>
typename ReactionNetwork<TImpl>::IndexType
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
    auto clusterData = _nw._clusterData;
    auto updater = typename Network::ClusterUpdater{};
    Kokkos::parallel_for(
            Range2D({0, 0}, {clusterData.numClusters, clusterData.gridSize}),
            KOKKOS_LAMBDA (IndexType i, IndexType j) {
        if (!xolotlCore::equal(clusterData.diffusionFactor(i), 0.0)) {
            updater.updateDiffusionCoefficient(clusterData, i, j);
        }
    });
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::generateClusterData(
    const typename Network::ClusterGenerator& generator)
{
    auto nClusters = _nw._clusterData.numClusters;
    _nw._clusterData.formationEnergy =
        Kokkos::View<double*>("Formation Energy", nClusters);
    _nw._clusterData.migrationEnergy =
        Kokkos::View<double*>("Migration Energy", nClusters);
    _nw._clusterData.diffusionFactor =
        Kokkos::View<double*>("Diffusion Factor", nClusters);
    _nw._clusterData.reactionRadius =
        Kokkos::View<double*>("Reaction Radius", nClusters);

    ClusterData data(_nw._clusterData);
    auto latticeParameter = _nw.getLatticeParameter();
    auto interstitialBias = _nw.getInterstitialBias();
    auto impurityRadius = _nw.getImpurityRadius();
    Kokkos::parallel_for(nClusters, KOKKOS_LAMBDA (const IndexType i) {
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
ReactionNetworkWorker<TImpl>::defineMomentIds()
{
    constexpr auto speciesRangeNoI = Network::getSpeciesRangeNoI();

    ClusterDataRef data(_nw._clusterData);

    auto nClusters = data.numClusters;
    auto counts = Kokkos::View<IndexType*>("Moment Id Counts", nClusters);

    IndexType nMomentIds = 0;
    Kokkos::parallel_reduce(nClusters,
        KOKKOS_LAMBDA (const IndexType i, IndexType& running)
        {
            const auto& reg = data.getCluster(i).getRegion();
            IndexType count = 0;
            for (auto k : speciesRangeNoI) {
                if (reg[k].length() != 1) {
                    ++count;
                }
            }
            running += count;
            counts(i) = count;
        },
        nMomentIds);

    Kokkos::parallel_scan(nClusters, KOKKOS_LAMBDA (IndexType i,
            IndexType& update, const bool finalPass) {
        const auto temp = counts(i);
        if (finalPass) {
            counts(i) = update;
        }
        update += temp;
    });

    Kokkos::parallel_for(nClusters, KOKKOS_LAMBDA (const IndexType i) {
        const auto& reg = data.getCluster(i).getRegion();
        IndexType current = counts(i);
        for (auto k : speciesRangeNoI) {
            if (reg[k].length() == 1) {
                data.momentIds(i, k()) = Network::invalidIndex();
            }
            else {
                data.momentIds(i, k()) = nClusters + current;
                ++current;
            }
        }
    });

    Kokkos::fence();
    _nw._numDOFs = nClusters + nMomentIds;
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::defineReactions()
{
    auto generator = _nw.asDerived()->getReactionGenerator();
    generator.generateReactions();

    _nw._reactionData = generator.getReactionData();
    _nw._reactions = ReactionCollection(generator.getProductionReactions(),
        generator.getDissociationReactions());
    _nw._sinkReactions = generator.getSinkReactions();
    _nw._resoReactions = generator.getReSoReactions();
}

template <typename TImpl>
typename ReactionNetworkWorker<TImpl>::IndexType
ReactionNetworkWorker<TImpl>::getDiagonalFill(
    typename Network::SparseFillMap& fillMap)
{
    auto connectivity = _nw._reactionData.connectivity;
    auto hConnRowMap = create_mirror_view(connectivity.row_map);
    deep_copy(hConnRowMap, connectivity.row_map);
    auto hConnEntries = create_mirror_view(connectivity.entries);
    deep_copy(hConnEntries, connectivity.entries);

    for (int i = 0; i < _nw.getDOF(); ++i) {
        auto jBegin = hConnRowMap(i);
        auto jEnd = hConnRowMap(i+1);
        std::vector<int> current;
        current.reserve(jEnd - jBegin);
        for (IndexType j = jBegin; j < jEnd; ++j) {
            current.push_back((int)hConnEntries(j));
        }
        fillMap[i] = std::move(current);
    }

    return hConnEntries.extent(0);
}
}
}
}
