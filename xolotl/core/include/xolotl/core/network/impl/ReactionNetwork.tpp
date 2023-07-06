#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/ReactionGenerator.h>
#include <xolotl/core/network/detail/TupleUtility.h>
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
	auto map = opts.getProcesses();
	this->setEnableStdReaction(map["reaction"]);
	this->setEnableReSolution(map["resolution"]);
	this->setEnableNucleation(map["heterogeneous"]);
	this->setEnableSink(map["sink"]);
	this->setEnableTrapMutation(map["modifiedTM"]);
	this->setEnableAttenuation(map["attenuation"]);
	this->setEnableConstantReaction(map["constant"]);
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

	// Skip the reactions for now if using constant reactions
	if (map["constant"])
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
ReactionNetwork<TImpl>::setEnableConstantReaction(bool reaction)
{
	this->_enableConstantReaction = reaction;
	_clusterData.h_view().setEnableConstantReaction(
		this->_enableConstantReaction);
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

	asDerived()->updateReactionRates(_currentTime);

	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setTime(double time)
{
	_currentTime = time;
	asDerived()->updateReactionRates(time);

	invalidateDataMirror();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::updateReactionRates(double time)
{
	_reactions.updateRates(time);
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
typename ReactionNetwork<TImpl>::MomentIdMap
ReactionNetwork<TImpl>::getAllMomentIdInfo()
{
	// Create the object to return
	MomentIdMap idMap;

	auto clusterDataMirror = this->getClusterDataMirror();

	// Loop on all the clusters
	for (auto i = 0; i < this->_numClusters; ++i) {
		auto cluster = clusterDataMirror.getCluster(i);
		auto momIds = cluster.getMomentIds();
		std::vector<IdType> temp;
		for (auto j = 0; j < momIds.extent(0); j++) {
			if (momIds(j) == this->invalidIndex())
				continue;
			temp.push_back(momIds(j));
		}
		idMap.push_back(temp);
	}
	return idMap;
}

template <typename TImpl>
std::string
ReactionNetwork<TImpl>::getHeaderString()
{
	// Create the object to return
	std::stringstream header;

	// Loop on all the clusters
	auto numSpecies = getSpeciesListSize();
	for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
		auto speciesName = this->getSpeciesName(id);
		header << speciesName << "_density " << speciesName << "_atom "
			   << speciesName << "_diameter " << speciesName
			   << "_partial_density " << speciesName << "_partial_atom "
			   << speciesName << "_partial_diameter ";
	}

	return header.str();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::initializeClusterMap(
	typename ReactionNetwork<TImpl>::BoundVector bounds,
	typename ReactionNetwork<TImpl>::MomentIdMapVector momIdInfos,
	typename ReactionNetwork<TImpl>::MomentIdMap fromSubNetwork)
{
	// Get the current bounds and moment id info
	auto currentBounds = getAllClusterBounds();
	auto currentMomIdInfo = getAllMomentIdInfo();

	// Check that the sizes add up
	IndexType nSubClusters = 0;
	for (auto subBounds : bounds) {
		nSubClusters += subBounds.size();
	}
	assert(this->_numClusters == nSubClusters);

	IndexType nSubMomIds = 0;
	for (auto subMomIdInfo : momIdInfos) {
		for (auto momIds : subMomIdInfo) {
			nSubMomIds += momIds.size();
		}
	}
	assert(this->_numDOFs == nSubClusters + nSubMomIds);

	// Create additional views to know how to treat each cluster
	auto dof = this->_numDOFs;

	for (auto subMap : fromSubNetwork) {
		// Get the sub DOF and initialize the rate map
		auto subDOF = subMap.size();
		auto hMap =
			Kokkos::View<IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
				subMap.data(), subDOF);
		auto dMap = Kokkos::View<IdType*>("Sub Map", subDOF);
		deep_copy(dMap, hMap);

		auto localIsInSub = BelongingView("In Sub", dof);
		auto localBackMap = OwnedSubMapView("Back Map", dof);

		// Initialize them
		Kokkos::parallel_for(
			dof, KOKKOS_LAMBDA(const IndexType i) {
				localIsInSub(i) = false;
				localBackMap(i) = 0;
			});
		Kokkos::parallel_for(
			dMap.extent(0), KOKKOS_LAMBDA(const IndexType i) {
				auto id = dMap(i);
				localIsInSub(id) = true;
				localBackMap(id) = i;
			});

		isInSub.push_back(localIsInSub);
		backMap.push_back(localBackMap);
	}

	return;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::initializeReactions()
{
	Connectivity connectivity;
	defineReactions(connectivity);
	generateDiagonalFill(connectivity);

	return;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setConstantRates(
	typename ReactionNetwork<TImpl>::RateVector rates)
{
	asDerived()->setConstantRates(rates);

	return;
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::setConstantConnectivities(
	typename ReactionNetwork<TImpl>::ConnectivitiesVector conns)
{
	asDerived()->setConstantConnectivities(conns);

	return;
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
void
ReactionNetwork<TImpl>::computeConstantRates(ConcentrationsView concentrations,
	RatesView rates, IndexType subId, IndexType gridIndex, double surfaceDepth,
	double spacing)
{
	asDerived()->computeConstantRatesPreProcess(
		concentrations, gridIndex, surfaceDepth, spacing);

	auto localInSub = isInSub[subId];
	auto localBackMap = backMap[subId];
	_reactions.forEach(DEVICE_LAMBDA(auto&& reaction) {
		reaction.contributeConstantRates(
			concentrations, rates, localInSub, localBackMap, gridIndex);
	});
	Kokkos::fence();
}

template <typename TImpl>
void
ReactionNetwork<TImpl>::getConstantConnectivities(
	ConnectivitiesView conns, IndexType subId)
{
	auto localInSub = isInSub[subId];
	auto localBackMap = backMap[subId];
	_reactions.forEach(DEVICE_LAMBDA(auto&& reaction) {
		reaction.contributeConstantConnectivities(
			conns, localInSub, localBackMap);
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

template <typename TReactionNetwork, typename TDerived>
struct TQMethodBase
{
	using AmountType = typename TReactionNetwork::AmountType;
	using Species = typename TReactionNetwork::Species;
	using Region = typename TReactionNetwork::Region;
	using ConcentrationsView = typename TReactionNetwork::ConcentrationsView;
	using ClusterData = typename TReactionNetwork::ClusterData;
	using TotalQuantity = typename TReactionNetwork::TotalQuantity;

	TQMethodBase(const TotalQuantity& quant) :
		_species(quant.speciesId.template cast<Species>()),
		_minSize(quant.minSize)
	{
	}

	KOKKOS_INLINE_FUNCTION
	void
	operator()(double concentration, double reactionRadius, const Region& clReg,
		double& lsum) const
	{
		using util::max;
		const auto& ival = clReg[_species];
		auto factor = clReg.volume() / ival.length();
		asDerived()->updateFactor(factor, concentration, reactionRadius);
		for (auto j = max(_minSize, ival.begin()); j < ival.end(); ++j) {
			lsum += asDerived()->contribution(j, factor);
		}
	}

	KOKKOS_INLINE_FUNCTION
	const TDerived*
	asDerived() const
	{
		return static_cast<const TDerived*>(this);
	}

	Species _species;
	AmountType _minSize;
};

template <typename TReactionNetwork>
struct TQMethodTotal :
	TQMethodBase<TReactionNetwork, TQMethodTotal<TReactionNetwork>>
{
	using Superclass =
		TQMethodBase<TReactionNetwork, TQMethodTotal<TReactionNetwork>>;
	using AmountType = typename Superclass::AmountType;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	void
	updateFactor(double& factor, double conc, double) const
	{
		factor *= conc;
	}

	KOKKOS_INLINE_FUNCTION
	double
	contribution(AmountType j, double factor) const
	{
		return factor;
	}
};

template <typename TReactionNetwork>
struct TQMethodAtom :
	TQMethodBase<TReactionNetwork, TQMethodAtom<TReactionNetwork>>
{
	using Superclass =
		TQMethodBase<TReactionNetwork, TQMethodAtom<TReactionNetwork>>;
	using AmountType = typename Superclass::AmountType;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	void
	updateFactor(double& factor, double conc, double) const
	{
		factor *= conc;
	}

	KOKKOS_INLINE_FUNCTION
	double
	contribution(AmountType j, double factor) const
	{
		return j * factor;
	}
};

template <typename TReactionNetwork>
struct TQMethodRadius :
	TQMethodBase<TReactionNetwork, TQMethodRadius<TReactionNetwork>>
{
	using Superclass =
		TQMethodBase<TReactionNetwork, TQMethodRadius<TReactionNetwork>>;
	using typename Superclass::AmountType;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	void
	updateFactor(double& factor, double conc, double reactionRadius) const
	{
		factor *= conc * reactionRadius;
	}

	KOKKOS_INLINE_FUNCTION
	double
	contribution(AmountType j, double factor) const
	{
		return factor;
	}
};

template <typename TReactionNetwork>
struct TQMethodVolume :
	TQMethodBase<TReactionNetwork, TQMethodVolume<TReactionNetwork>>
{
	using Superclass =
		TQMethodBase<TReactionNetwork, TQMethodVolume<TReactionNetwork>>;
	using typename Superclass::AmountType;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	void
	updateFactor(double& factor, double conc, double reactionRadius) const
	{
		factor *= conc * reactionRadius * reactionRadius * reactionRadius;
	}

	KOKKOS_INLINE_FUNCTION
	double
	contribution(AmountType j, double factor) const
	{
		return factor;
	}
};

template <typename TReactionNetwork>
struct TQMethodTrapped :
	TQMethodBase<TReactionNetwork, TQMethodTrapped<TReactionNetwork>>
{
	using Superclass =
		TQMethodBase<TReactionNetwork, TQMethodTrapped<TReactionNetwork>>;
	using typename Superclass::AmountType;
	using typename Superclass::Region;
	using typename Superclass::Species;

	using Superclass::Superclass;

	struct VacancyDetect
	{
		bool hasVacancy;
		Species vIndex;
	};

	static constexpr VacancyDetect
	checkVacancy() noexcept
	{
		for (auto s : TReactionNetwork::getSpeciesRangeNoI()) {
			if (isVacancy(s)) {
				return {true, s};
			}
		}
		return {false, Species{}};
	}

	static constexpr VacancyDetect vd = checkVacancy();

	KOKKOS_INLINE_FUNCTION
	void
	operator()(double concentration, double reactionRadius, const Region& clReg,
		double& lsum) const
	{
		if constexpr (vd.hasVacancy) {
			if (clReg[vd.vIndex].begin() > 0) {
				Superclass::operator()(
					concentration, reactionRadius, clReg, lsum);
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	updateFactor(double& factor, double conc, double reactionRadius) const
	{
		factor *= conc;
	}

	KOKKOS_INLINE_FUNCTION
	double
	contribution(AmountType j, double factor) const
	{
		return j * factor;
	}
};

template <typename... Ts>
struct TQMethodChain
{
	template <typename T>
	TQMethodChain(T&&)
	{
	}

	template <typename... TArgs>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(TArgs&&...) const
	{
	}
};

template <typename THead, typename... TTail>
struct TQMethodChain<THead, TTail...> : TQMethodChain<TTail...>
{
	using Superclass = TQMethodChain<TTail...>;
	using Method = THead;
	using TotalQuantity = typename Method::TotalQuantity;
	using ClusterData = typename Method::ClusterData;
	using Region = typename Method::Region;

	static constexpr std::size_t _index = sizeof...(TTail);

	TQMethodChain(const TotalQuantity quantities[]) :
		Superclass(quantities),
		_method(quantities[_index])
	{
	}

	KOKKOS_INLINE_FUNCTION
	void
	operator()(double concentration, double reactionRadius, const Region& clReg,
		double dst[]) const
	{
		Superclass::operator()(concentration, reactionRadius, clReg, dst);

		_method(concentration, reactionRadius, clReg, dst[_index]);
	}

	Method _method;
};

template <typename TReactionNetwork, typename... TQMethods>
class TotalQuantityReduceFunctor
{
public:
	static constexpr std::size_t N = sizeof...(TQMethods);

	// Required by Kokkos
	using execution_space = Kokkos::DefaultExecutionSpace;
	using value_type = double[];
	using size_type = typename execution_space::size_type;
	const unsigned value_count = N;

	using ConcentrationsView = typename TReactionNetwork::ConcentrationsView;
	using TilesView = typename TReactionNetwork::Subpaving::TilesView;
	using ClusterDataView = typename TReactionNetwork::ClusterDataView;
	using TotalQuantity = typename TReactionNetwork::TotalQuantity;

	TotalQuantityReduceFunctor(ConcentrationsView concs, TilesView tiles,
		ClusterDataView clData,
		const util::Array<TotalQuantity, N>& quantities) :
		_concentrations(concs),
		_tiles(tiles),
		_clusterData(clData),
		_methods(quantities.data())
	{
	}

	KOKKOS_INLINE_FUNCTION
	void
	init(double dst[]) const
	{
		for (unsigned i = 0; i < value_count; ++i) {
			dst[i] = 0.0;
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	join(double dst[], const double src[]) const
	{
		for (unsigned i = 0; i < value_count; ++i) {
			dst[i] += src[i];
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	operator()(size_type i, double dst[]) const
	{
		_methods(_concentrations[i], _clusterData().reactionRadius[i],
			_tiles[i].getRegion(), dst);
	}

private:
	ConcentrationsView _concentrations;
	TilesView _tiles;
	ClusterDataView _clusterData;

	using ReverseMethodsTuple = detail::TupleReverse<std::tuple<TQMethods...>>;
	using MethodChain =
		detail::TupleApplyAll<TQMethodChain, ReverseMethodsTuple>;
	MethodChain _methods;
};

template <std::size_t N>
struct TQArrayMatcher
{
	using TotalQuantity = IReactionNetwork::TotalQuantity;
	using Q = TotalQuantity::Type;

	const util::Array<TotalQuantity, N>& quantities;

	bool
	operator()(const util::Array<Q, N>& check) const
	{
		for (std::size_t i = 0; i < N; ++i) {
			if (quantities[i].type != check[i]) {
				return false;
			}
		}
		return true;
	}
};

template <typename TImpl>
template <typename... TQMethods>
util::Array<double, sizeof...(TQMethods)>
ReactionNetwork<TImpl>::getTotalsImpl(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, sizeof...(TQMethods)>& quantities)
{
	constexpr std::size_t N = sizeof...(TQMethods);

	auto tiles = _subpaving.getTiles();
	auto clusterData = _clusterData.d_view;
	auto result = util::Array<double, N>{};

	auto temp =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
			result.data(), result.size());

	Kokkos::parallel_reduce("ReactionNetworkWorker::getTotals",
		this->_numClusters,
		TotalQuantityReduceFunctor<TImpl, TQMethods...>{
			concentrations, tiles, clusterData, quantities},
		temp);

	Kokkos::fence();

	return result;
}

struct TotalQuantityError : std::runtime_error
{
	using TotalQuantity = IReactionNetwork::TotalQuantity;

	explicit TotalQuantityError(const std::string& message) :
		std::runtime_error(
			"Unsupported total quantity specification ( " + message + " )")
	{
	}

	template <typename TArray>
	TotalQuantityError(const std::string& message, TArray&& quantities) :
		TotalQuantityError(message + "\ngiven: " + toString(quantities))
	{
	}

	static std::string
	toString(TotalQuantity::Type type)
	{
		static constexpr const char* qTypes[] = {
			"total", "atom", "radius", "volume", "trapped"};
		using Id = std::underlying_type_t<TotalQuantity::Type>;
		return qTypes[static_cast<Id>(type)];
	}

	template <typename TArray>
	static std::string
	toString(TArray&& quantities)
	{
		std::ostringstream oss;
		oss << "[" << toString(quantities[0].type);
		for (std::size_t i = 1; i < quantities.size(); ++i) {
			oss << ", " << toString(quantities[i].type);
		}
		oss << "]";
		return oss.str();
	}
};

template <typename TImpl>
util::Array<double, 1>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 1>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<1>{quantities};
	if (match({{Q::total}})) {
		return getTotalsImpl<TQMethodTotal>(concentrations, quantities);
	}
	if (match({{Q::atom}})) {
		return getTotalsImpl<TQMethodAtom>(concentrations, quantities);
	}

	throw TotalQuantityError(
		"getTotals<1>: must specify one of { [total], [atom] }", quantities);
}

template <typename TImpl>
util::Array<double, 2>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 2>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<2>{quantities};
	if (match({{Q::atom, Q::atom}})) {
		return getTotalsImpl<TQMethodAtom, TQMethodAtom>(
			concentrations, quantities);
	}

	throw TotalQuantityError(
		"getTotals<2>: must specify [atom, atom]", quantities);
}

template <typename TImpl>
util::Array<double, 3>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 3>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<3>{quantities};
	if (match({{Q::atom, Q::atom, Q::atom}})) {
		return getTotalsImpl<TQMethodAtom, TQMethodAtom, TQMethodAtom>(
			concentrations, quantities);
	}

	throw TotalQuantityError(
		"getTotals<3>: must specify [atom, atom, atom]", quantities);
}

template <typename TImpl>
util::Array<double, 4>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 4>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<4>{quantities};
	if (match({{Q::total, Q::radius, Q::total, Q::radius}})) {
		return getTotalsImpl<TQMethodTotal, TQMethodRadius, TQMethodTotal,
			TQMethodRadius>(concentrations, quantities);
	}
	if (match({{Q::atom, Q::atom, Q::atom, Q::atom}})) {
		return getTotalsImpl<TQMethodAtom, TQMethodAtom, TQMethodAtom,
			TQMethodAtom>(concentrations, quantities);
	}

	throw TotalQuantityError("getTotals<4>: must specify one of {\n"
							 "\t[total, radius, total, radius],\n"
							 "\t[atom, atom, atom, atom]\n"
							 "}",
		quantities);
}

template <typename TImpl>
util::Array<double, 5>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 5>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<5>{quantities};
	if (match({{Q::atom, Q::atom, Q::atom, Q::atom, Q::atom}})) {
		return getTotalsImpl<TQMethodAtom, TQMethodAtom, TQMethodAtom,
			TQMethodAtom, TQMethodAtom>(concentrations, quantities);
	}

	throw TotalQuantityError(
		"getTotals<5>: must specify [atom, atom, atom, atom, atom]",
		quantities);
}

template <typename TImpl>
util::Array<double, 6>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 6>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<6>{quantities};
	if (match({{Q::total, Q::atom, Q::radius, Q::total, Q::atom, Q::radius}})) {
		return getTotalsImpl<TQMethodTotal, TQMethodAtom, TQMethodRadius,
			TQMethodTotal, TQMethodAtom, TQMethodRadius>(
			concentrations, quantities);
	}
	if (match(
			{{Q::total, Q::atom, Q::radius, Q::total, Q::radius, Q::volume}})) {
		return getTotalsImpl<TQMethodTotal, TQMethodAtom, TQMethodRadius,
			TQMethodTotal, TQMethodRadius, TQMethodVolume>(
			concentrations, quantities);
	}

	throw TotalQuantityError("getTotals<6>: must specify one of{\n"
							 "\t[total, atom, radius, total, atom, radius],\n"
							 "\t[total, atom, radius, total, radius, volume]\n"
							 "}",
		quantities);
}

template <typename TImpl>
util::Array<double, 7>
ReactionNetwork<TImpl>::getTotals(ConcentrationsView concentrations,
	const util::Array<TotalQuantity, 7>& quantities)
{
	using Q = TotalQuantity::Type;
	auto match = TQArrayMatcher<7>{quantities};
	if (match({{Q::total, Q::atom, Q::radius, Q::total, Q::atom, Q::radius,
			Q::volume}})) {
		return getTotalsImpl<TQMethodTotal, TQMethodAtom, TQMethodRadius,
			TQMethodTotal, TQMethodAtom, TQMethodRadius, TQMethodVolume>(
			concentrations, quantities);
	}

	throw TotalQuantityError(
		"getTotals<7>: must specify\n"
		"\t[total, atom, radius, total, atom, radius, volume]",
		quantities);
}

template <typename TImpl>
std::vector<double>
ReactionNetwork<TImpl>::getTotalsVec(ConcentrationsView concentrations,
	const std::vector<TotalQuantity>& quantities)
{
	auto toArray = [&quantities](auto&& a) {
		std::copy(begin(quantities), end(quantities), begin(a));
		return a;
	};
	auto toVector = [](auto&& a) {
		std::vector<double> result;
		result.reserve(a.size());
		std::copy(begin(a), end(a), std::back_inserter(result));
		return result;
	};
	switch (quantities.size()) {
	case 1:
		using A1 = util::Array<TotalQuantity, 1>;
		return toVector(getTotals(concentrations, toArray(A1{})));
	case 2:
		using A2 = util::Array<TotalQuantity, 2>;
		return toVector(getTotals(concentrations, toArray(A2{})));
	case 3:
		using A3 = util::Array<TotalQuantity, 3>;
		return toVector(getTotals(concentrations, toArray(A3{})));
	case 4:
		using A4 = util::Array<TotalQuantity, 4>;
		return toVector(getTotals(concentrations, toArray(A4{})));
	case 5:
		using A5 = util::Array<TotalQuantity, 5>;
		return toVector(getTotals(concentrations, toArray(A5{})));
	case 6:
		using A6 = util::Array<TotalQuantity, 6>;
		return toVector(getTotals(concentrations, toArray(A6{})));
	default:
		throw TotalQuantityError(
			"getTotalsVec: Currently supports only up to 6; given size = " +
			std::to_string(quantities.size()));
		break;
	}
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
	generator.setConstantConnectivities(_nw._constantConns);
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
