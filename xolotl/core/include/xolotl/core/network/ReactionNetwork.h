#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <type_traits>

#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <xolotl/core/network/Cluster.h>
#include <xolotl/core/network/IReactionNetwork.h>
#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>
#include <xolotl/core/network/detail/ReactionCollection.h>
#include <xolotl/options/IOptions.h>
#include <xolotl/util/NotImplementedError.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TImpl, typename TDerived>
class ReactionGeneratorBase;
} // namespace detail

template <typename TImpl>
struct ReactionNetworkInterface
{
	using Type = IReactionNetwork;
};

template <typename TImpl>
class ReactionNetwork : public ReactionNetworkInterface<TImpl>::Type
{
	template <typename, typename>
	friend class detail::ReactionGeneratorBase;

public:
	using Superclass = typename ReactionNetworkInterface<TImpl>::Type;
	using Traits = ReactionNetworkTraits<TImpl>;
	using Species = typename Traits::Species;

private:
	static_assert(std::is_base_of<IReactionNetwork, Superclass>::value,
		"ReactionNetwork must inherit from IReactionNetwork");

	using Types = detail::ReactionNetworkTypes<TImpl>;

	static constexpr std::size_t numSpecies = Traits::numSpecies;

public:
	using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
	using SpeciesRange = EnumSequenceRange<Species, numSpecies>;
	using ClusterGenerator = typename Traits::ClusterGenerator;
	using ClusterUpdater = typename Types::ClusterUpdater;
	using Connectivity = typename Superclass::Connectivity;
	using AmountType = typename IReactionNetwork::AmountType;
	using IndexType = typename IReactionNetwork::IndexType;
	using Subpaving = typename Types::Subpaving;
	using SubpavingMirror = typename Subpaving::HostMirror;
	using SubdivisionRatio = plsm::SubdivisionRatio<numSpecies>;
	using Composition = typename Types::Composition;
	using Region = typename Types::Region;
	using Ival = typename Region::IntervalType;
	using ConcentrationsView = typename IReactionNetwork::ConcentrationsView;
	using FluxesView = typename IReactionNetwork::FluxesView;
	using RatesView = typename IReactionNetwork::RatesView;
	using ConnectivitiesView = typename IReactionNetwork::ConnectivitiesView;
	using ConnectivitiesPairView =
		typename IReactionNetwork::ConnectivitiesPairView;
	using SubMapView = typename IReactionNetwork::SubMapView;
	using OwnedSubMapView = typename IReactionNetwork::OwnedSubMapView;
	using BelongingView = typename IReactionNetwork::BelongingView;
	using SparseFillMap = typename IReactionNetwork::SparseFillMap;
	using ClusterData = typename Types::ClusterData;
	using ClusterDataMirror = typename Types::ClusterDataMirror;
	using ClusterDataView = Kokkos::View<ClusterData>;
	using ClusterDataHostView = typename ClusterDataView::host_mirror_type;
	using ReactionCollection = typename Types::ReactionCollection;
	using Bounds = IReactionNetwork::Bounds;
	using BoundVector = IReactionNetwork::BoundVector;
	using MomentIdMap = IReactionNetwork::MomentIdMap;
	using MomentIdMapVector = IReactionNetwork::MomentIdMapVector;
	using RateVector = IReactionNetwork::RateVector;
	using ConnectivitiesPair = IReactionNetwork::ConnectivitiesPair;
	using PhaseSpace = IReactionNetwork::PhaseSpace;
	using TotalQuantity = IReactionNetwork::TotalQuantity;

	template <typename PlsmContext>
	using Cluster = Cluster<TImpl, PlsmContext>;

	void
	copyClusterDataView();

	ReactionNetwork() = default;

	ReactionNetwork(const Subpaving& subpaving, IndexType gridSize,
		const options::IOptions& opts);

	ReactionNetwork(const Subpaving& subpaving, IndexType gridSize);

	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		const std::vector<SubdivisionRatio>& subdivisionRatios,
		IndexType gridSize, const options::IOptions& opts);

	ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		IndexType gridSize, const options::IOptions& opts);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::size_t
	getNumberOfSpecies() noexcept
	{
		return SpeciesSequence::size();
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr std::size_t
	getNumberOfSpeciesNoI() noexcept
	{
		return SpeciesSequence::sizeNoI();
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr SpeciesRange
	getSpeciesRange() noexcept
	{
		return SpeciesRange{};
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr SpeciesRange
	getSpeciesRangeNoI() noexcept
	{
		return SpeciesRange(
			SpeciesSequence::first(), SpeciesSequence::lastNoI());
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr SpeciesRange
	getSpeciesRangeForGrouping() noexcept
	{
		using GroupingRange =
			SpeciesForGrouping<Species, SpeciesSequence::size()>;
		return SpeciesRange(GroupingRange::first, GroupingRange::last);
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<Species>
	mapToMomentId(EnumSequence<Species, SpeciesSequence::size()> value)
	{
		using GroupingRange =
			SpeciesForGrouping<Species, SpeciesSequence::size()>;
		return GroupingRange::mapToMomentId(value);
	}

	void
	initializeExtraClusterData(const options::IOptions&)
	{
	}

	void
	updateExtraClusterData(
		const std::vector<double>&, const std::vector<double>&)
	{
	}

	std::size_t
	getSpeciesListSize() const noexcept override
	{
		return getNumberOfSpecies();
	}

	const std::string&
	getSpeciesLabel(SpeciesId id) const override;

	const std::string&
	getSpeciesName(SpeciesId id) const override;

	SpeciesId
	parseSpeciesId(const std::string& speciesLabel) const override;

	void
	setLatticeParameter(double latticeParameter) override;

	void
	setImpurityRadius(double impurityRadius) noexcept override
	{
		this->_impurityRadius =
			asDerived()->checkImpurityRadius(impurityRadius);
	}

	void
	setFissionRate(double rate) override;

	void
	setZeta(double zeta) override;

	void
	setEnableStdReaction(bool reaction) override;

	void
	setEnableReSolution(bool reaction) override;

	void
	setEnableNucleation(bool reaction) override;

	void
	setEnableSink(bool reaction) override;

	void
	setEnableTrapMutation(bool reaction) override;

	void
	setEnableConstantReaction(bool reaction) override;

	void
	setEnableReducedJacobian(bool reduced) override;

	void
	setEnableReadRates(bool read) override;

	void
	setGridSize(IndexType gridSize) override;

	void
	setTemperatures(const std::vector<double>& gridTemperatures,
		const std::vector<double>& gridDepths) override;

	void
	setTime(double time) override;

	const std::vector<AmountType>&
	getMinRadiusSizes() const override
	{
		return _minRadiusSizes;
	}

	[[noreturn]] std::string
	getMonitorOutputFileName() const override
	{
		throw util::NotImplementedError();
	}

	[[noreturn]] std::string
	getMonitorDataHeaderString() const override
	{
		throw util::NotImplementedError();
	}

	std::vector<double>
	getMonitorDataValues(Kokkos::View<const double*> conc, double fac) override
	{
		auto ret = std::vector<double>(getMonitorDataLineSize(), 0.0);
		addMonitorDataValues(conc, fac, ret);
		return ret;
	}

	void
	addMonitorDataValues(Kokkos::View<const double*> conc, double fac,
		std::vector<double>& totalVals) override
	{
		throw util::NotImplementedError();
	}

	[[noreturn]] std::size_t
	getMonitorDataLineSize() const override
	{
		throw util::NotImplementedError();
	}

	void
	writeMonitorOutputHeader() const override
	{
		std::ofstream(this->getMonitorOutputFileName())
			<< this->getMonitorDataHeaderString() << std::endl;
	}

	void
	writeMonitorDataLine(
		const std::vector<double>& localData, double time) override
	{
		throw util::NotImplementedError();
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept override;

	void
	syncClusterDataOnHost() override;

	void
	invalidateDataMirror()
	{
		_subpavingMirror.reset();
		_clusterDataMirror.reset();
	}

	const SubpavingMirror&
	getSubpavingMirror()
	{
		if (!_subpavingMirror.has_value()) {
			syncClusterDataOnHost();
		}
		return *_subpavingMirror;
	}

	const ClusterDataMirror&
	getClusterDataMirror()
	{
		if (!_clusterDataMirror.has_value()) {
			syncClusterDataOnHost();
		}
		return *_clusterDataMirror;
	}

	template <typename MemSpace>
	KOKKOS_INLINE_FUNCTION
	Cluster<MemSpace>
	findCluster(const Composition& comp, MemSpace)
	{
		if constexpr (!std::is_same_v<plsm::HostMemSpace,
						  plsm::DeviceMemSpace> &&
			std::is_same_v<MemSpace, plsm::HostMemSpace>) {
#ifdef __CUDACC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20011, 20014
#endif
			// CUDA warns about calling __host__ functions here.
			// However, his branch is never compiled using CUDA
			auto id = getSubpavingMirror().findTileId(comp);
			return getClusterDataMirror().getCluster(
				id == _subpaving.invalidIndex() ? this->invalidIndex() :
												  IndexType(id));
#ifdef __CUDACC__
#pragma nv_diagnostic pop
#endif
		}
		else {
			auto id = _subpaving.findTileId(comp);
			return _clusterData.d_view().getCluster(
				id == _subpaving.invalidIndex() ? this->invalidIndex() :
												  IndexType(id));
		}
	}

	KOKKOS_INLINE_FUNCTION
	auto
	findCluster(const Composition& comp)
	{
		return findCluster(comp, plsm::DeviceMemSpace{});
	}

	IndexType
	findClusterId(const std::vector<AmountType>& composition) override
	{
		assert(composition.size() == getNumberOfSpecies());
		Composition comp;
		for (std::size_t i = 0; i < composition.size(); ++i) {
			comp[i] = composition[i];
		}
		return findCluster(comp, plsm::HostMemSpace{}).getId();
	}

	ClusterCommon<plsm::DeviceMemSpace>
	getClusterCommon(IndexType clusterId, plsm::DeviceMemSpace) override
	{
		return ClusterCommon<plsm::DeviceMemSpace>(
			_clusterData.d_view.data(), clusterId);
	}

	ClusterCommon<plsm::HostMemSpace>
	getClusterCommon(IndexType clusterId) override
	{
		return getClusterDataMirror().getClusterCommon(clusterId);
	}

	IndexType
	getLargestClusterId() override
	{
		return asDerived()->checkLargestClusterId();
	}

	Bounds
	getAllClusterBounds() override;

	MomentIdMap
	getAllMomentIdInfo() override;

	void initializeClusterMap(
		BoundVector, MomentIdMapVector, MomentIdMap) override;

	void
	initializeReactions() override;

	void
	setConstantRates(RatesView, IndexType gridIndex) override;

	void setConstantConnectivities(ConnectivitiesPair) override;

	void
	initializeRateEntries(const ConnectivitiesPair&, IndexType) override;

	void
	initializeRateEntries(const std::vector<ConnectivitiesPair>&) override;

	void
	setConstantRateEntries() override;

	PhaseSpace
	getPhaseSpace() override;

	template <typename MemSpace>
	KOKKOS_INLINE_FUNCTION
	Cluster<MemSpace>
	getCluster(IndexType clusterId, MemSpace)
	{
		if constexpr (!std::is_same_v<plsm::HostMemSpace,
						  plsm::DeviceMemSpace> &&
			std::is_same_v<MemSpace, plsm::HostMemSpace>) {
#ifdef __CUDACC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20011, 20014
#endif
			// CUDA warns about calling a __host__ function here.
			// However, his branch is never compiled using CUDA
			return getClusterDataMirror().getCluster(clusterId);
#ifdef __CUDACC__
#pragma nv_diagnostic pop
#endif
		}
		else {
			return _clusterData.d_view().getCluster(clusterId);
		}
	}

	KOKKOS_INLINE_FUNCTION
	auto
	getCluster(IndexType clusterId)
	{
		return getCluster(clusterId, plsm::DeviceMemSpace{});
	}

	KOKKOS_INLINE_FUNCTION
	Subpaving&
	getSubpaving()
	{
		return _subpaving;
	}

	void
	computeFluxesPreProcess(
		ConcentrationsView, FluxesView, IndexType, double, double)
	{
	}

	void
	computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex = 0, double surfaceDepth = 0.0,
		double spacing = 0.0) final;

	template <typename TReaction>
	void
	computeFluxes(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex = 0, double surfaceDepth = 0.0,
		double spacing = 0.0)
	{
		asDerived()->computeFluxesPreProcess(
			concentrations, fluxes, gridIndex, surfaceDepth, spacing);

		_reactions.template forEachOn<TReaction>(
			"ReactionNetwork::computeFluxes", DEVICE_LAMBDA(auto&& reaction) {
				reaction.contributeFlux(concentrations, fluxes, gridIndex);
			});
		Kokkos::fence();
	}

	void
	computePartialsPreProcess(
		ConcentrationsView, FluxesView, IndexType, double, double)
	{
	}

	void
	computeAllPartials(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex = 0,
		double surfaceDepth = 0.0, double spacing = 0.0) override;

	void
	computeConstantRatesPreProcess(
		ConcentrationsView, IndexType, double, double)
	{
	}

	void
	computeConstantRates(ConcentrationsView concentrations, RatesView rates,
		IndexType subId, IndexType gridIndex = 0, double surfaceDepth = 0.0,
		double spacing = 0.0) final;

	void
	getConstantConnectivities(ConnectivitiesView conns, IndexType subId) final;

	template <typename TReaction>
	void
	computePartials(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex = 0,
		double surfaceDepth = 0.0, double spacing = 0.0)
	{
		// Reset the values
		const auto& nValues = values.extent(0);
		Kokkos::parallel_for(
			"ReactionNetwork::computePartials::resetValues", nValues,
			KOKKOS_LAMBDA(const IndexType i) { values(i) = 0.0; });

		asDerived()->computePartialsPreProcess(
			concentrations, values, gridIndex, surfaceDepth, spacing);

		if (this->_enableReducedJacobian) {
			_reactions.template forEachOn<TReaction>(
				"ReactionNetwork::computePartials",
				DEVICE_LAMBDA(auto&& reaction) {
					reaction.contributeReducedPartialDerivatives(
						concentrations, values, gridIndex);
				});
		}
		else {
			_reactions.template forEachOn<TReaction>(
				"ReactionNetwork::computePartials",
				DEVICE_LAMBDA(auto&& reaction) {
					reaction.contributePartialDerivatives(
						concentrations, values, gridIndex);
				});
		}

		Kokkos::fence();
	}

	double
	getLargestRate() override;

	double
	getLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex) override;

	IndexType
	getDiagonalFill(SparseFillMap& fillMap) override;

	template <typename... TQMethods>
	util::Array<double, sizeof...(TQMethods)>
	getTotalsImpl(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, sizeof...(TQMethods)>& quantities);

	template <template <typename> typename... TQMethods>
	util::Array<double, sizeof...(TQMethods)>
	getTotalsImpl(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, sizeof...(TQMethods)>& quantities)
	{
		return getTotalsImpl<TQMethods<TImpl>...>(concentrations, quantities);
	}

	util::Array<double, 1>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 1>& quantities) override;

	util::Array<double, 2>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 2>& quantities) override;

	util::Array<double, 3>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 3>& quantities) override;

	util::Array<double, 4>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 4>& quantities) override;

	util::Array<double, 5>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 5>& quantities) override;

	util::Array<double, 6>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 6>& quantities) override;

	util::Array<double, 7>
	getTotals(ConcentrationsView concentrations,
		const util::Array<TotalQuantity, 7>& quantities) override;

	std::vector<double>
	getTotalsVec(ConcentrationsView concentrations,
		const std::vector<TotalQuantity>& quantities) override;

	/**
	 * Get the total concentration of a given type of clusters.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @return The total accumulated concentration
	 */
	double
	getTotalConcentration(ConcentrationsView concentrations, Species type,
		AmountType minSize = 0);

	double
	getTotalConcentration(ConcentrationsView concentrations, SpeciesId species,
		AmountType minSize = 0) override
	{
		auto type = species.cast<Species>();
		return getTotalConcentration(concentrations, type, minSize);
	}

	/**
	 * Get the total concentration of a given type of clusters times their
	 * radius.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @return The total accumulated concentration times the radius of the
	 * cluster
	 */
	double
	getTotalRadiusConcentration(ConcentrationsView concentrations, Species type,
		AmountType minSize = 0);

	double
	getTotalRadiusConcentration(ConcentrationsView concentrations,
		SpeciesId species, AmountType minSize = 0) override
	{
		auto type = species.cast<Species>();
		return getTotalRadiusConcentration(concentrations, type, minSize);
	}

	/**
	 * Get the total concentration of a given type of clusters times the number
	 * of atoms.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @return The total accumulated concentration times the size of the cluster
	 */
	double
	getTotalAtomConcentration(ConcentrationsView concentrations, Species type,
		AmountType minSize = 0);

	double
	getTotalAtomConcentration(ConcentrationsView concentrations,
		SpeciesId species, AmountType minSize = 0) override
	{
		auto type = species.cast<Species>();
		return getTotalAtomConcentration(concentrations, type, minSize);
	}

	/**
	 * Get the averaged species to defect ratio.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @return The ratio
	 */
	double
	getTotalVolumeRatio(ConcentrationsView concentrations, Species type,
		AmountType minSize = 0);

	double
	getTotalVolumeRatio(ConcentrationsView concentrations, SpeciesId species,
		AmountType minSize = 0) override
	{
		auto type = species.cast<Species>();
		return getTotalVolumeRatio(concentrations, type, minSize);
	}

	/**
	 * Get the variance associated with averaged species to defect ratio.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @param mean The ration mean
	 * @return The variance
	 */
	double
	getTotalRatioVariance(ConcentrationsView concentrations, Species type,
		double mean, AmountType minSize = 0);

	double
	getTotalRatioVariance(ConcentrationsView concentrations, SpeciesId species,
		double mean, AmountType minSize = 0) override
	{
		auto type = species.cast<Species>();
		return getTotalRatioVariance(concentrations, type, minSize);
	}

	/**
	 * Get the total concentration of a given type of clusters only if it is
	 * trapped in a vacancy.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @return The total accumulated concentration times the size of the cluster
	 */
	double
	getTotalTrappedAtomConcentration(ConcentrationsView concentrations,
		Species type, AmountType minSize = 0);

	/**
	 * Get the total concentration of a given type of clusters times their
	 * volume.
	 *
	 * @param concentration The vector of concentrations
	 * @param type The type of atom we want the concentration of
	 * @param minSize The minimum number of atom to start counting
	 * @return The total accumulated volume fraction
	 */
	double
	getTotalVolumeFraction(ConcentrationsView concentrations, Species type,
		AmountType minSize = 0);

	void
	updateReactionRates(double time = 0.0);

	void
	updateOutgoingDiffFluxes(double* gridPointSolution, double factor,
		std::vector<IndexType> diffusingIds, std::vector<double>& fluxes,
		IndexType gridIndex) override;

	void
	updateOutgoingAdvecFluxes(double* gridPointSolution, double factor,
		std::vector<IndexType> advectingIds, std::vector<double> sinkStrengths,
		std::vector<double>& fluxes, IndexType gridIndex) override;

private:
	KOKKOS_INLINE_FUNCTION
	TImpl*
	asDerived()
	{
		return static_cast<TImpl*>(this);
	}

	static std::map<std::string, SpeciesId>
	createSpeciesLabelMap() noexcept;

	static std::vector<AmountType>
	computeMinRadiusSizes(const options::IOptions& opts);

	void
	defineMomentIds();

	void
	generateClusterData(const ClusterGenerator& generator);

	void
	readClusters(const std::string filename = "reactionRates.txt");

	void
	readReactions(
		double temperature, const std::string filename = "reactionRates.txt");

	void
	defineReactions(Connectivity& connectivity);

	void
	updateDiffusionCoefficients();

	KOKKOS_INLINE_FUNCTION
	double
	getTemperature(IndexType gridIndex) const noexcept
	{
		return _clusterData.d_view().temperature(gridIndex);
	}

private:
	void
	generateDiagonalFill(const Connectivity& connectivity);

private:
	std::optional<SubpavingMirror> _subpavingMirror;

	SparseFillMap _connectivityMap;

	std::vector<BelongingView> isInSub;
	std::vector<OwnedSubMapView> backMap;

protected:
	std::optional<ClusterDataMirror> _clusterDataMirror;

	Kokkos::DualView<ClusterData> _clusterData;

	Subpaving _subpaving;

	ReactionCollection _reactions;

	std::map<std::string, SpeciesId> _speciesLabelMap;

	// Reaction energies
	Kokkos::View<double**> _reactionEnergies;

	ConnectivitiesPairView _constantConnsRows;
	ConnectivitiesPairView _constantConnsEntries;

	double _currentTime;

	std::vector<AmountType> _minRadiusSizes;
};

namespace detail
{
template <typename TImpl>
class DefaultClusterUpdater
{
public:
	using Network = ReactionNetwork<TImpl>;
	using Types = ReactionNetworkTypes<TImpl>;
	using ClusterData = typename Types::ClusterData;
	using IndexType = typename Types::IndexType;

	KOKKOS_INLINE_FUNCTION
	void
	updateDiffusionCoefficient(const ClusterData& data, IndexType clusterId,
		IndexType gridIndex) const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
