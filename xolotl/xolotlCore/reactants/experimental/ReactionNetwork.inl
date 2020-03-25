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

    this->_atomicVolume =
        0.5 * lParam * lParam * lParam;
    auto mirror = Kokkos::create_mirror_view(_clusterData.atomicVolume);
    mirror(0) = this->_atomicVolume;
    Kokkos::deep_copy(_clusterData.atomicVolume, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setGridSize(IndexType gridSize)
{
    _gridSize = gridSize;
    _clusterData.setGridSize(_gridSize);
    _clusterDataMirror.setGridSize(_gridSize);
    _reactionData.setGridSize(_gridSize);

    auto prodReactions = _prodReactions;
    auto dissReactions = _dissReactions;
    auto numProdReactions = _numProdReactions;
    auto reactionData = _reactionData;
    auto clusterData = _clusterData;
    Kokkos::parallel_for(_numReactions, KOKKOS_LAMBDA (IndexType i) {
        if (i < numProdReactions) {
            prodReactions(i).updateData(reactionData, clusterData, i);
        }
        else {
            dissReactions(i-numProdReactions).updateData(reactionData, clusterData, i);
        }
    });
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setTemperatures(const std::vector<double>& gridTemps)
{
    Kokkos::View<const double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
        tempsHost(gridTemps.data(), this->_gridSize);
    Kokkos::deep_copy(_clusterData.temperature, tempsHost);

    updateDiffusionCoefficients();

    auto prodReactions = _prodReactions;
    auto dissReactions = _dissReactions;
    auto numProdReactions = _numProdReactions;
    Kokkos::parallel_for(_numReactions, KOKKOS_LAMBDA (IndexType i) {
        if (i < numProdReactions) {
            prodReactions(i).updateRates();
        }
        else {
            dissReactions(i-numProdReactions).updateRates();
        }
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
    // Get the extent of the reactions view
	auto prodReactions = _prodReactions;
    auto dissReactions = _dissReactions;
    auto numProdReactions = _numProdReactions;
    auto numReactions = _numReactions;
    // Loop on the reactions
    Kokkos::parallel_for(numReactions, KOKKOS_LAMBDA (const IndexType i) {
        if (i < numProdReactions) {
            prodReactions(i).contributeFlux(concentrations, fluxes, gridIndex);
        }
        else {
            dissReactions(i-numProdReactions).contributeFlux(concentrations, fluxes, gridIndex);
        }
    });
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

    // Get the extent of the reactions view
    auto prodReactions = _prodReactions;
    auto dissReactions = _dissReactions;
    auto numProdReactions = _numProdReactions;
    auto connectivity = _reactionData.connectivity;
    // Loop on the reactions
    Kokkos::parallel_for(_numReactions, KOKKOS_LAMBDA (const IndexType i) {
        if (i < numProdReactions) {
            prodReactions(i).contributePartialDerivatives(concentrations,
                values, connectivity, gridIndex);
        }
        else {
            dissReactions(i-numProdReactions).contributePartialDerivatives(concentrations,
                values, connectivity, gridIndex);
        }
    });
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
    auto prodReactions = _prodReactions;
    auto dissReactions = _dissReactions;
    auto numProdReactions = _numProdReactions;
    double leftSideRate = 0.0;
    // Loop on all the rates to get the maximum
    Kokkos::parallel_reduce(_numReactions,
            KOKKOS_LAMBDA (const IndexType i, double &lsum) {
        if (i < numProdReactions) {
            lsum += prodReactions(i).contributeLeftSideRate(concentrations,
                clusterId, gridIndex);
        }
        else {
            lsum += dissReactions(i-numProdReactions).contributeLeftSideRate(concentrations,
                clusterId, gridIndex);
        }
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
    auto data = _nw._clusterData;
    Kokkos::parallel_for(Range2D({0, 0}, {data.numClusters, data.gridSize}),
            KOKKOS_LAMBDA (IndexType i, IndexType j) {
        if (!xolotlCore::equal(data.diffusionFactor(i), 0.0))
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
    _nw._prodReactions = generator.getProductionReactions();
    _nw._dissReactions = generator.getDissociationReactions();
    _nw._numProdReactions = _nw._prodReactions.size();
    _nw._numReactions = _nw._numProdReactions + _nw._dissReactions.size();
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
