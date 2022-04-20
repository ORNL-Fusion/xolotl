#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/ReactionGenerator.h>
#include <xolotl/core/network/detail/impl/ClusterData.tpp>
#include <xolotl/core/network/detail/impl/ReactionGenerator.tpp>
#include <xolotl/core/network/impl/Reaction.tpp>
#include <xolotl/options/Options.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TImpl>
inline void
ReactionNetwork<TImpl>::copyClusterDataView()
{
	_clusterData.modify_host();
	_clusterData.sync_device();
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(const Subpaving& subpaving,
	IndexType gridSize, const options::IOptions& opts) :
	Superclass(gridSize),
	_subpaving(subpaving),
	_clusterData("Cluster Data"),
	_worker(*this),
	_speciesLabelMap(createSpeciesLabelMap())
{
	_clusterData.h_view() = ClusterData(_subpaving, gridSize);
	copyClusterDataView();

	this->setMaterial(opts.getMaterial());

	// Set constants
	this->setInterstitialBias(opts.getBiasFactor());
	this->setImpurityRadius(opts.getImpurityRadius());
	this->setLatticeParameter(opts.getLatticeParameter());
	this->setFissionRate(opts.getFluxAmplitude());
	this->setZeta(opts.getZeta());
	this->setTauBursting(opts.getBurstingDepth());
	this->setFBursting(opts.getBurstingFactor());
	auto map = opts.getProcesses();
	this->setEnableStdReaction(map["reaction"]);
	this->setEnableReSolution(map["resolution"]);
	this->setEnableNucleation(map["heterogeneous"]);
	this->setEnableSink(map["sink"]);
	this->setEnableTrapMutation(map["modifiedTM"]);
	this->setEnableAttenuation(map["attenuation"]);
	this->setEnableBursting(map["bursting"]);
	std::string petscString = opts.getPetscArg();
	auto tokens = util::Tokenizer<>{petscString}();
	bool useReduced = false;
	for (const auto& token : tokens) {
		if (token == "-snes_mf_operator") {
			useReduced = true;
			break;
		}
	}
	this->setEnableReducedJacobian(useReduced);

	this->_numClusters = _clusterData.h_view().numClusters;
	asDerived()->initializeExtraClusterData(opts);
	generateClusterData(ClusterGenerator{opts});
	defineMomentIds();

	if (map["largeBubble"])
		return;

	Connectivity connectivity;
	defineReactions(connectivity);
	generateDiagonalFill(connectivity);
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
	const Subpaving& subpaving, IndexType gridSize) :
	ReactionNetwork(subpaving, gridSize, options::Options{})
{
}

template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts) :
	ReactionNetwork(
		[&]() -> Subpaving {
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
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts) :
	ReactionNetwork(
		maxSpeciesAmounts,
		[&maxSpeciesAmounts]() -> std::vector<SubdivisionRatio> {
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
const std::string&
ReactionNetwork<TImpl>::getSpeciesLabel(SpeciesId id) const
{
	return toLabelString(id.cast<Species>());
}

template <typename TImpl>
const std::string&
ReactionNetwork<TImpl>::getSpeciesName(SpeciesId id) const
{
	return toNameString(id.cast<Species>());
}

template <typename TImpl>
SpeciesId
ReactionNetwork<TImpl>::parseSpeciesId(const std::string& speciesLabel) const
{
	auto it = _speciesLabelMap.find(speciesLabel);
	if (it == _speciesLabelMap.end()) {
		throw InvalidSpeciesId("Unrecognized species type: " + speciesLabel);
	}
	return it->second;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setLatticeParameter(double latticeParameter)
{
	auto lParam = asDerived()->checkLatticeParameter(latticeParameter);
	this->_latticeParameter = lParam;
	_clusterData.h_view().setLatticeParameter(this->_latticeParameter);

	this->_atomicVolume = asDerived()->computeAtomicVolume(lParam);
	_clusterData.h_view().setAtomicVolume(this->_atomicVolume);

	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setFissionRate(double rate)
{
	Superclass::setFissionRate(rate);
	_clusterData.h_view().setFissionRate(this->_fissionRate);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setZeta(double z)
{
	_clusterData.h_view().setZeta(z);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setTauBursting(double tau)
{
	_clusterData.h_view().setTauBursting(tau);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setFBursting(double f)
{
	_clusterData.h_view().setFBursting(f);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableStdReaction(bool reaction)
{
	Superclass::setEnableStdReaction(reaction);
	_clusterData.h_view().setEnableStdReaction(this->_enableStdReaction);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableReSolution(bool reaction)
{
	Superclass::setEnableReSolution(reaction);
	_clusterData.h_view().setEnableReSolution(this->_enableReSolution);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableNucleation(bool reaction)
{
	Superclass::setEnableNucleation(reaction);
	_clusterData.h_view().setEnableNucleation(this->_enableNucleation);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableSink(bool reaction)
{
	this->_enableSink = reaction;
	_clusterData.h_view().setEnableSink(this->_enableSink);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableTrapMutation(bool reaction)
{
	Superclass::setEnableTrapMutation(reaction);
	_clusterData.h_view().setEnableTrapMutation(this->_enableTrapMutation);
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableBursting(bool reaction)
{
	this->_enableBursting = reaction;
	_clusterData.h_view().setEnableBurst(this->_enableBursting);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setEnableReducedJacobian(bool reduced)
{
	this->_enableReducedJacobian = reduced;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setGridSize(IndexType gridSize)
{
	this->_gridSize = gridSize;
	// Set value in mirror so we don't need to invalidate
	if (_clusterDataMirror.has_value()) {
		_clusterDataMirror.value().setGridSize(gridSize);
	}
	_clusterData.h_view().setGridSize(gridSize);
	copyClusterDataView();
	_reactions.setGridSize(gridSize);
	_reactions.updateAll(_clusterData.d_view);
	Kokkos::fence();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setTemperatures(
	const std::vector<double>& gridTemps, const std::vector<double>& gridDepths)
{
	Kokkos::View<const double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		tempsHost(gridTemps.data(), this->_gridSize);
	Kokkos::deep_copy(_clusterData.h_view().temperature, tempsHost);

	updateDiffusionCoefficients();

	asDerived()->updateExtraClusterData(gridTemps, gridDepths);

	asDerived()->updateReactionRates();

	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::updateReactionRates()
{
	_reactions.updateRates();
}

template <typename TImpl>
std::uint64_t
ReactionNetwork<TImpl>::getDeviceMemorySize() const noexcept
{
	std::uint64_t ret = _subpaving.getDeviceMemorySize();

	ret += _clusterData.h_view().getDeviceMemorySize();
	ret += _reactions.getDeviceMemorySize();

	return ret;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::syncClusterDataOnHost()
{
	BOOST_LOG_FUNCTION();
	XOLOTL_LOG_XTRA;

	_subpavingMirror = _subpaving.makeMirrorCopy();

	auto dataMirror = ClusterDataMirror(*_subpavingMirror, this->_gridSize);
	dataMirror.deepCopy(_clusterData.h_view());
	_clusterDataMirror = dataMirror;
}

template <typename TImpl>
ClusterCommon<plsm::HostMemSpace>
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
	if (hasVacancy)
		comp[vIndex] = 1;

	auto clusterId = findCluster(comp, plsm::HostMemSpace{}).getId();

	return getClusterDataMirror().getClusterCommon(clusterId);
}

template <typename TImpl>
typename ReactionNetwork<TImpl>::Bounds
ReactionNetwork<TImpl>::getAllClusterBounds()
{
	// Create the object to return
	Bounds bounds;

	// Loop on all the clusters
	constexpr auto speciesRange = getSpeciesRange();
	auto tiles = getSubpavingMirror().getTiles();
	for (IndexType i = 0; i < this->_numClusters; ++i) {
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
		space.push_back(toLabelString(j));
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
	_clusterData.h_view().generate(generator, this->getLatticeParameter(),
		this->getInterstitialBias(), this->getImpurityRadius());
	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::computeAllFluxes(ConcentrationsView concentrations,
	FluxesView fluxes, IndexType gridIndex, double surfaceDepth, double spacing)
{
	asDerived()->computeFluxesPreProcess(
		concentrations, fluxes, gridIndex, surfaceDepth, spacing);

	_reactions.forEach(
		"ReactionNetwork::computeAllFluxes", DEVICE_LAMBDA(auto&& reaction) {
			reaction.contributeFlux(concentrations, fluxes, gridIndex);
		});
	Kokkos::fence();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::computeAllPartials(ConcentrationsView concentrations,
	Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
	double spacing)
{
	// Reset the values
	const auto& nValues = values.extent(0);
	Kokkos::parallel_for(
		"ReactionNetwork::computeAllPartials::resetValues", nValues,
		KOKKOS_LAMBDA(const IndexType i) { values(i) = 0.0; });

	asDerived()->computePartialsPreProcess(
		concentrations, values, gridIndex, surfaceDepth, spacing);

	if (this->_enableReducedJacobian) {
		_reactions.forEach(
			"ReactionNetwork::computeAllPartials",
			DEVICE_LAMBDA(auto&& reaction) {
				reaction.contributeReducedPartialDerivatives(
					concentrations, values, gridIndex);
			});
	}
	else {
		_reactions.forEach(
			"ReactionNetwork::computeAllPartials",
			DEVICE_LAMBDA(auto&& reaction) {
				reaction.contributePartialDerivatives(
					concentrations, values, gridIndex);
			});
	}

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
ReactionNetwork<TImpl>::getLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// Get the extent of the reactions
	double leftSideRate = 0.0;
	// Loop on all the rates to get the maximum
	_reactions.reduce(
		"ReactionNetwork::getLeftSideRate",
		DEVICE_LAMBDA(auto&& reaction, double& lsum) {
			lsum += reaction.contributeLeftSideRate(
				concentrations, clusterId, gridIndex);
		},
		leftSideRate);

	return leftSideRate;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize)
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
ReactionNetwork<TImpl>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize)
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
	if (!hasVacancy)
		return 0.0;

	auto tiles = _subpaving.getTiles();
	double conc = 0.0;
	Kokkos::parallel_reduce(
		"ReactionNetwork::getTotalTrappedAtomConcentration", this->_numClusters,
		KOKKOS_LAMBDA(IndexType i, double& lsum) {
			const Region& clReg = tiles(i).getRegion();
			if (clReg[vIndex].begin() > 0) {
				const auto factor = clReg.volume() / clReg[type].length();
				for (AmountType j : makeIntervalRange(clReg[type])) {
					if (j >= minSize)
						lsum += concentrations(i) * j * factor;
				}
			}
		},
		conc);

	Kokkos::fence();

	return conc;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::updateOutgoingDiffFluxes(double* gridPointSolution,
	double factor, std::vector<IndexType> diffusingIds,
	std::vector<double>& fluxes, IndexType gridIndex)
{
	// Loop on the diffusing clusters
	for (auto l : diffusingIds) {
		// Get the cluster and composition
		auto cluster = this->getClusterCommon(l);
		auto reg = this->getCluster(l, plsm::HostMemSpace{}).getRegion();
		Composition comp = reg.getOrigin();
		// Get its concentration
		double conc = gridPointSolution[l];
		// Get its size and diffusion coefficient
		int size = 0;
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			auto type = id.cast<Species>();
			size += comp[type];
		}
		double coef = cluster.getDiffusionCoefficient(gridIndex);
		// Compute the flux
		double newFlux = (double)size * factor * coef * conc;

		// Check the cluster type
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			auto type = id.cast<Species>();
			if (comp.isOnAxis(type)) {
				fluxes[id()] += newFlux;
				break;
			}
		}
	}

	return;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::updateOutgoingAdvecFluxes(double* gridPointSolution,
	double factor, std::vector<IndexType> advectingIds,
	std::vector<double> sinkStrengths, std::vector<double>& fluxes,
	IndexType gridIndex)
{
	int advClusterIdx = 0;
	// Loop on the advecting clusters
	for (auto l : advectingIds) {
		// Get the cluster and composition
		auto cluster = this->getClusterCommon(l);
		auto reg = this->getCluster(l, plsm::HostMemSpace{}).getRegion();
		Composition comp = reg.getOrigin();
		// Get its concentration
		double conc = gridPointSolution[l];
		// Get its size and diffusion coefficient
		int size = 0;
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			auto type = id.cast<Species>();
			size += comp[type];
		}
		double coef = cluster.getDiffusionCoefficient(gridIndex);
		// Compute the flux if the temperature is valid
		double newFlux = 0.0;
		if (cluster.getTemperature(gridIndex) > 0.0)
			newFlux = (double)size * factor * coef * conc *
				sinkStrengths[advClusterIdx] /
				cluster.getTemperature(gridIndex);

		// Check the cluster type
		for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
			auto type = id.cast<Species>();
			if (comp.isOnAxis(type)) {
				fluxes[id()] += newFlux;
				break;
			}
		}

		advClusterIdx++;
	}

	return;
}

template <typename TImpl>
double
ReactionNetwork<TImpl>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize)
{
	return _worker.getTotalVolumeFraction(concentrations, type, minSize);
}

template <typename TImpl>
std::map<std::string, SpeciesId>
ReactionNetwork<TImpl>::createSpeciesLabelMap() noexcept
{
	std::map<std::string, SpeciesId> labelMap;
	for (auto s : getSpeciesRange()) {
		labelMap.emplace(
			toLabelString(s.value), SpeciesId(s.value, getNumberOfSpecies()));
	}
	return labelMap;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineMomentIds()
{
	_worker.defineMomentIds();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::defineReactions(Connectivity& connectivity)
{
	_worker.defineReactions(connectivity);
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::generateDiagonalFill(const Connectivity& connectivity)
{
	auto hConnRowMap = create_mirror_view(connectivity.row_map);
	deep_copy(hConnRowMap, connectivity.row_map);
	auto hConnEntries = create_mirror_view(connectivity.entries);
	deep_copy(hConnEntries, connectivity.entries);

	_connectivityMap.clear();
	for (int i = 0; i < this->getDOF(); ++i) {
		auto jBegin = hConnRowMap(i);
		auto jEnd = hConnRowMap(i + 1);
		std::vector<int> current;
		current.reserve(jEnd - jBegin);
		for (IndexType j = jBegin; j < jEnd; ++j) {
			current.push_back((int)hConnEntries(j));
		}
		_connectivityMap[i] = std::move(current);
	}
}

template <typename TImpl>
typename ReactionNetwork<TImpl>::IndexType
ReactionNetwork<TImpl>::getDiagonalFill(SparseFillMap& fillMap)
{
	IndexType nnz = 0;
	for (int i = 0; i < this->getDOF(); ++i) {
		const auto& current = _connectivityMap[i];
		nnz += current.size();
		fillMap.insert_or_assign(i, current);
	}
	return nnz;
}

namespace detail
{
template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::updateDiffusionCoefficients()
{
	using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
	auto clusterData = _nw._clusterData.d_view;
	auto updater = typename Network::ClusterUpdater{};
	Kokkos::parallel_for(
		"ReactionNetworkWorker::updateDiffusionCoefficients",
		Range2D({0, 0},
			{_nw._clusterData.h_view().numClusters,
				_nw._clusterData.h_view().gridSize}),
		KOKKOS_LAMBDA(IndexType i, IndexType j) {
			if (!util::equal(clusterData().diffusionFactor(i), 0.0)) {
				updater.updateDiffusionCoefficient(clusterData(), i, j);
			}
		});
	Kokkos::fence();
	_nw.invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::defineMomentIds()
{
	constexpr auto speciesRange = Network::getSpeciesRangeForGrouping();

	auto nClusters = _nw._clusterData.h_view().numClusters;
	auto counts = Kokkos::View<IndexType*>("Moment Id Counts", nClusters);

	auto data = _nw._clusterData.d_view.data();

	IndexType nMomentIds = 0;
	Kokkos::parallel_reduce(
		"ReactionNetworkWorker::defineMomentIds::count", nClusters,
		KOKKOS_LAMBDA(const IndexType i, IndexType& running) {
			const auto& reg = data->getCluster(i).getRegion();
			IndexType count = 0;
			for (auto k : speciesRange) {
				if (reg[k].length() != 1) {
					++count;
				}
			}
			running += count;
			counts(i) = count;
		},
		nMomentIds);

	Kokkos::parallel_scan(
		"ReactionNetworkWorker::defineMomentIds::scan", nClusters,
		KOKKOS_LAMBDA(IndexType i, IndexType & update, const bool finalPass) {
			const auto temp = counts(i);
			if (finalPass) {
				counts(i) = update;
			}
			update += temp;
		});

	Kokkos::parallel_for(
		"ReactionNetworkWorker::defineMomentIds::assignMomentIds", nClusters,
		KOKKOS_LAMBDA(const IndexType i) {
			for (IndexType j = 0; j < data->momentIds.extent(1); ++j) {
				data->momentIds(i, j) = Network::invalidIndex();
			}

			const auto& reg = data->getCluster(i).getRegion();
			IndexType current = counts(i);
			for (auto k : speciesRange) {
				if (reg[k].length() == 1) {
					if (data->momentIds(i, Network::mapToMomentId(k)) ==
						nClusters + current - 1)
						continue;
					data->momentIds(i, Network::mapToMomentId(k)) =
						Network::invalidIndex();
				}
				else {
					data->momentIds(i, Network::mapToMomentId(k)) =
						nClusters + current;
					++current;
				}
			}
		});

	Kokkos::fence();
	_nw._numDOFs = nClusters + nMomentIds;
	_nw.invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetworkWorker<TImpl>::defineReactions(Connectivity& connectivity)
{
	auto generator = _nw.asDerived()->getReactionGenerator();
	_nw._reactions = generator.generateReactions();
	connectivity = generator.getConnectivity();
}

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize)
{
	auto tiles = _nw._subpaving.getTiles();
	double conc = 0.0;
	Kokkos::parallel_reduce(
		"ReactionNetworkWorker::getTotalConcentration", _nw._numClusters,
		KOKKOS_LAMBDA(IndexType i, double& lsum) {
			using util::max;
			const auto& clReg = tiles(i).getRegion();
			const auto& ival = clReg[type];
			const auto factor = clReg.volume() / ival.length();
			for (auto j = max(minSize, ival.begin()); j < ival.end(); ++j) {
				lsum += concentrations(i) * factor;
			}
		},
		conc);

	Kokkos::fence();

	return conc;
}

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize)
{
	auto tiles = _nw._subpaving.getTiles();
	double conc = 0.0;
	auto clusterData = _nw._clusterData.d_view;
	Kokkos::parallel_reduce(
		"ReactionNetworkWorker::getTotalRadiusConcentration", _nw._numClusters,
		KOKKOS_LAMBDA(IndexType i, double& lsum) {
			using util::max;
			const auto& clReg = tiles(i).getRegion();
			const auto& ival = clReg[type];
			const auto factor = clReg.volume() / ival.length();
			for (auto j = max(minSize, ival.begin()); j < ival.end(); ++j) {
				lsum += concentrations(i) * clusterData().reactionRadius(i) *
					factor;
			}
		},
		conc);

	Kokkos::fence();

	return conc;
}

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize)
{
	auto tiles = _nw._subpaving.getTiles();
	double conc = 0.0;
	Kokkos::parallel_reduce(
		"ReactionNetworkWorker::getTotalAtomConcentration", _nw._numClusters,
		KOKKOS_LAMBDA(IndexType i, double& lsum) {
			using util::max;
			const auto& clReg = tiles(i).getRegion();
			const auto& ival = clReg[type];
			const auto factor = clReg.volume() / ival.length();
			for (auto j = max(minSize, ival.begin()); j < ival.end(); ++j) {
				lsum += concentrations(i) * j * factor;
			}
		},
		conc);

	Kokkos::fence();

	return conc;
}

template <typename TImpl>
double
ReactionNetworkWorker<TImpl>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize)
{
	auto tiles = _nw._subpaving.getTiles();
	double conc = 0.0;
	auto clusterData = _nw._clusterData.d_view;
	Kokkos::parallel_reduce(
		"ReactionNetworkWorker::getTotalVolumeFraction", _nw._numClusters,
		KOKKOS_LAMBDA(IndexType i, double& lsum) {
			using util::max;
			const auto& clReg = tiles(i).getRegion();
			const auto& ival = clReg[type];
			const auto factor = clReg.volume() / ival.length();
			const auto rRad = clusterData().reactionRadius(i);
			const auto rRad3 = rRad * rRad * rRad;
			for (auto j = max(minSize, ival.begin()); j < ival.end(); ++j) {
				lsum += concentrations(i) * rRad3 * factor;
			}
		},
		conc);

	Kokkos::fence();

	constexpr double sphereFactor = 4.0 * ::xolotl::core::pi / 3.0;

	return conc * sphereFactor;
}

template <typename TImpl>
KOKKOS_INLINE_FUNCTION
void
DefaultClusterUpdater<TImpl>::updateDiffusionCoefficient(
	const ClusterData& data, IndexType clusterId, IndexType gridIndex) const
{
	data.diffusionCoefficient(clusterId, gridIndex) =
		data.diffusionFactor(clusterId) *
		exp(-data.migrationEnergy(clusterId) /
			(kBoltzmann * data.temperature(gridIndex)));
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
