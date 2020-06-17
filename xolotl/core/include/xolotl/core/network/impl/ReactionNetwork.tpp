#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/ReactionGenerator.h>
#include <xolotl/core/network/detail/impl/ReactionGenerator.tpp>
#include <xolotl/core/network/impl/Reaction.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(const Subpaving& subpaving,
        IndexType gridSize, const options::IOptions& opts)
    :
    IReactionNetwork(gridSize),
    _subpaving(subpaving),
    _clusterData(_subpaving, gridSize),
    _worker(*this)
{
    // Set constants
    setInterstitialBias(opts.getBiasFactor());
    setImpurityRadius(opts.getImpurityRadius());
    setLatticeParameter(opts.getLatticeParameter());
    setFissionRate(opts.getFluxAmplitude());
    setZeta(opts.getZeta());
    auto map = opts.getProcesses();
    setEnableStdReaction(map["reaction"]);
    setEnableReSolution(map["resolution"]);
    setEnableNucleation(map["heterogeneous"]);

    auto tiles = subpaving.getTiles(plsm::onDevice);
    _numClusters = tiles.extent(0);
    
    asDerived()->checkTiles(opts);
    
    // PRINT ALL THE CLUSTERS
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
//    std::cout << "num: " << _numClusters << std::endl;

    generateClusterData(ClusterGenerator{opts});
    defineMomentIds();

    defineReactions();
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(const Subpaving& subpaving,
        IndexType gridSize)
    :
    ReactionNetwork(subpaving, gridSize, options::Options{})
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
        const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        IndexType gridSize, const options::IOptions& opts)
    :
    ReactionNetwork(
        [&]() -> Subpaving
        {
            Region latticeRegion{};
            for (std::size_t i = 0; i < getNumberOfSpecies(); ++i) {
                latticeRegion[i] = Ival{0, maxSpeciesAmounts[i] + 1};
            }
            Subpaving sp(latticeRegion, subdivisionRatios);
            sp.refine(ClusterGenerator{opts});
            return sp;
        }(),
        gridSize, opts)
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
        const std::vector<AmountType>& maxSpeciesAmounts,
        IndexType gridSize, const options::IOptions& opts)
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
        gridSize, opts)
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

    this->_atomicVolume = asDerived()->computeAtomicVolume(lParam);
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
ReactionNetwork<TImpl>::setZeta(double z)
{
    auto mirror = Kokkos::create_mirror_view(_clusterData.zeta);
    mirror(0) = z;
    Kokkos::deep_copy(_clusterData.zeta, mirror);
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
ReactionNetwork<TImpl>::setEnableReSolution(bool reaction)
{
    this->_enableReSolution = reaction;
    auto mirror = Kokkos::create_mirror_view(_clusterData.enableReSolution);
    mirror(0) = this->_enableReSolution;
    Kokkos::deep_copy(_clusterData.enableReSolution, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableNucleation(bool reaction)
{
    this->_enableNucleation = reaction;
    auto mirror = Kokkos::create_mirror_view(_clusterData.enableNucleation);
    mirror(0) = this->_enableNucleation;
    Kokkos::deep_copy(_clusterData.enableNucleation, mirror);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setGridSize(IndexType gridSize)
{
    _gridSize = gridSize;
    _clusterData.setGridSize(_gridSize);
    _clusterDataMirror.setGridSize(_gridSize);
    _reactions.setGridSize(_gridSize);
    _reactions.updateAll(_clusterData);
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
    Kokkos::fence();
}

template <typename TImpl>
std::uint64_t
ReactionNetwork<TImpl>::getDeviceMemorySize() const noexcept
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

    ret += _reactions.getDeviceMemorySize();

    return ret;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::syncClusterDataOnHost()
{
    _subpaving.syncTiles(plsm::onHost);
    auto mirror = ClusterDataMirror(_subpaving, this->_gridSize);
    Kokkos::deep_copy(mirror.atomicVolume, _clusterData.atomicVolume);
    Kokkos::deep_copy(mirror.latticeParameter, _clusterData.latticeParameter);
    Kokkos::deep_copy(mirror.fissionRate, _clusterData.fissionRate);
    Kokkos::deep_copy(mirror.enableStdReaction, _clusterData.enableStdReaction);
    Kokkos::deep_copy(mirror.enableReSolution, _clusterData.enableReSolution);
    Kokkos::deep_copy(mirror.enableNucleation, _clusterData.enableNucleation);
    Kokkos::deep_copy(mirror.temperature, _clusterData.temperature);
    Kokkos::deep_copy(mirror.momentIds, _clusterData.momentIds);
    Kokkos::deep_copy(mirror.reactionRadius, _clusterData.reactionRadius);
    Kokkos::deep_copy(mirror.formationEnergy, _clusterData.formationEnergy);
    Kokkos::deep_copy(mirror.migrationEnergy, _clusterData.migrationEnergy);
    Kokkos::deep_copy(mirror.diffusionFactor, _clusterData.diffusionFactor);
    Kokkos::deep_copy(mirror.diffusionCoefficient, _clusterData.diffusionCoefficient);
    _clusterDataMirror = mirror;
}

template <typename TImpl>
KOKKOS_INLINE_FUNCTION
typename ReactionNetwork<TImpl>::template Cluster<plsm::OnDevice>
ReactionNetwork<TImpl>::findCluster(const Composition& comp,
    plsm::OnDevice context)
{
    return Cluster<plsm::OnDevice>(_clusterData,
        _subpaving.findTileId(comp, context));
}

template <typename TImpl>
typename ReactionNetwork<TImpl>::template Cluster<plsm::OnHost>
ReactionNetwork<TImpl>::findCluster(const Composition& comp,
    plsm::OnHost context)
{
    return Cluster<plsm::OnHost>(_clusterDataMirror,
        _subpaving.findTileId(comp, context));
}

template <typename TImpl>
ClusterCommon<plsm::OnHost>
ReactionNetwork<TImpl>::getSingleVacancy()
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

template <typename TImpl>
typename ReactionNetwork<TImpl>::Bounds
ReactionNetwork<TImpl>::getAllClusterBounds()
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

template <typename TImpl>
typename ReactionNetwork<TImpl>::PhaseSpace
ReactionNetwork<TImpl>::getPhaseSpace()
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

    auto connectivity = _reactions.getConnectivity();
    _reactions.apply(DEVICE_LAMBDA (auto&& reaction) {
        reaction.contributePartialDerivatives(concentrations, values,
            connectivity, gridIndex);
    });

    Kokkos::fence();
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getLargestRate()
{
    return _reactions.getLargestRate();
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
    return _worker.getTotalConcentration(concentrations, type, minSize);
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalRadiusConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize)
{
    return _worker.getTotalRadiusConcentration(concentrations, type, minSize);
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize)
{
    return _worker.getTotalAtomConcentration(concentrations, type, minSize);
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

    Kokkos::fence();

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
        if (!util::equal(clusterData.diffusionFactor(i), 0.0)) {
            updater.updateDiffusionCoefficient(clusterData, i, j);
        }
    });
    Kokkos::fence();
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
        data.diffusionFactor(i) = generator.getDiffusionFactor(cluster, latticeParameter);
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
    _nw._reactions = generator.generateReactions();

    // _nw._reactionData = generator.getReactionData();
    // _nw._reactions = generator.getReactionCollection();
}

template <typename TImpl>
typename ReactionNetworkWorker<TImpl>::IndexType
ReactionNetworkWorker<TImpl>::getDiagonalFill(
    typename Network::SparseFillMap& fillMap)
{
    auto connectivity = _nw._reactions.getConnectivity();
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

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize)
{
    auto tiles = _nw._subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    Kokkos::parallel_reduce(_nw._numClusters,
            KOKKOS_LAMBDA (IndexType i, double& lsum) {
    	const auto& clReg = tiles(i).getRegion();
    	const auto factor = clReg.volume() / clReg[type].length();
    	for (AmountType j : makeIntervalRange(clReg[type])) {
    		if (j >= minSize) lsum += concentrations(i) * factor;
    	}
    }, conc);

    Kokkos::fence();

    return conc;
}

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalRadiusConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize)
{
    auto tiles = _nw._subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    auto clusterData = _nw._clusterData;
    Kokkos::parallel_reduce(_nw._numClusters,
            KOKKOS_LAMBDA (IndexType i, double& lsum) {
    	const auto& clReg = tiles(i).getRegion();
    	const auto factor = clReg.volume() / clReg[type].length();
    	for (AmountType j : makeIntervalRange(clReg[type])) {
    		if (j >= minSize) 
    			lsum += concentrations(i) * clusterData.reactionRadius(i) * factor;
    	}
    }, conc);

    Kokkos::fence();

    return conc;
}

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalAtomConcentration(
    ConcentrationsView concentrations, Species type, AmountType minSize)
{
    auto tiles = _nw._subpaving.getTiles(plsm::onDevice);
    double conc = 0.0;
    Kokkos::parallel_reduce(_nw._numClusters,
            KOKKOS_LAMBDA (IndexType i, double &lsum) {
    	const auto& clReg = tiles(i).getRegion();
    	const auto factor = clReg.volume() / clReg[type].length();
    	for (AmountType j : makeIntervalRange(clReg[type])) {
    		if (j >= minSize) lsum += concentrations(i) * j * factor;
    	}
    }, conc);

    Kokkos::fence();

    return conc;
}

template <typename TImpl>
KOKKOS_INLINE_FUNCTION
void
DefaultClusterUpdater<TImpl>::updateDiffusionCoefficient(
    const ClusterData& data, IndexType clusterId, IndexType gridIndex) const
{
    data.diffusionCoefficient(clusterId, gridIndex) =
        data.diffusionFactor(clusterId) * exp(
            -data.migrationEnergy(clusterId) /
            (kBoltzmann * data.temperature(gridIndex)));
}
}
}
}
}
