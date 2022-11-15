#pragma once

#include <string>
#include <unordered_map>

#include <Kokkos_Core.hpp>

#include <xolotl/core/network/Cluster.h>
#include <xolotl/core/network/SpeciesId.h>
#include <xolotl/core/network/detail/ClusterConnectivity.h>
#include <xolotl/core/network/detail/ClusterData.h>
#include <xolotl/core/network/detail/ReactionData.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Virtual class interface for accessing anything needed by outside
 * classes like the solver.
 */
class IReactionNetwork
{
public:
	using IndexType = detail::ReactionNetworkIndexType;
	using AmountType = detail::CompositionAmountType;
	using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
	using OwnedConcentrationsView = Kokkos::View<double*>;
	using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
	using OwnedFluxesView = Kokkos::View<double*>;
	using RatesView = Kokkos::View<double**>;
	using ConnectivitiesView = Kokkos::View<bool**>;
	using SubMapView = Kokkos::View<AmountType*, Kokkos::MemoryUnmanaged>;
	using OwnedSubMapView = Kokkos::View<AmountType*>;
	using BelongingView = Kokkos::View<bool*>;
	using Connectivity = detail::ClusterConnectivity<>;
	using SparseFillMap = std::unordered_map<int, std::vector<int>>;
	using Bounds = std::vector<std::vector<AmountType>>;
	using BoundVector = std::vector<std::vector<std::vector<AmountType>>>;
	using MomentIdMap = std::vector<std::vector<IdType>>;
	using MomentIdMapVector = std::vector<std::vector<std::vector<IdType>>>;
	using RateVector = std::vector<std::vector<double>>;
	using ConnectivitiesVector = std::vector<std::vector<bool>>;
	using PhaseSpace = std::vector<std::string>;

	KOKKOS_INLINE_FUNCTION
	static constexpr IndexType
	invalidIndex() noexcept
	{
		return detail::invalidNetworkIndex;
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr AmountType
	invalidAmount() noexcept
	{
		return detail::invalidSpeciesAmount;
	}

	IReactionNetwork() = default;

	IReactionNetwork(IndexType gridSize) : _gridSize(gridSize)
	{
	}

	virtual ~IReactionNetwork()
	{
	}

	virtual std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return 0;
	}

	virtual std::size_t
	getSpeciesListSize() const noexcept = 0;

	virtual SpeciesId
	parseSpeciesId(const std::string& label) const = 0;

	virtual const std::string&
	getSpeciesLabel(SpeciesId id) const = 0;

	virtual const std::string&
	getSpeciesName(SpeciesId id) const = 0;

	/**
	 * @brief Degrees of freedom includes _numClusters
	 * plus the moments.
	 */
	KOKKOS_INLINE_FUNCTION
	IndexType
	getDOF() const noexcept
	{
		return _numDOFs;
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	getNumClusters() const noexcept
	{
		return _numClusters;
	}

	const std::string&
	getMaterial() const noexcept
	{
		return _material;
	}

	void
	setMaterial(const std::string& mat)
	{
		_material = mat;
	}

	KOKKOS_INLINE_FUNCTION
	double
	getLatticeParameter() const noexcept
	{
		return _latticeParameter;
	}

	virtual void
	setLatticeParameter(double latticeParameter) = 0;

	KOKKOS_INLINE_FUNCTION
	double
	getAtomicVolume() const noexcept
	{
		return _atomicVolume;
	}

	KOKKOS_INLINE_FUNCTION
	double
	getInterstitialBias() const noexcept
	{
		return _interstitialBias;
	}

	void
	setInterstitialBias(double interstitialBias) noexcept
	{
		_interstitialBias = interstitialBias;
	}

	KOKKOS_INLINE_FUNCTION
	double
	getImpurityRadius() const noexcept
	{
		return _impurityRadius;
	}

	virtual void
	setImpurityRadius(double impurityRadius) noexcept
	{
		_impurityRadius = impurityRadius;
	}

	KOKKOS_INLINE_FUNCTION
	double
	getFissionRate() const noexcept
	{
		return _fissionRate;
	}

	virtual void
	setFissionRate(double rate)
	{
		_fissionRate = rate;
	}

	/**
	 * @brief Zeta is used to compute the rates of re-solution reactions.
	 */
	virtual void
	setZeta(double z) = 0;

	bool
	getEnableStdReaction() const noexcept
	{
		return _enableStdReaction;
	}

	virtual void
	setEnableStdReaction(bool reaction)
	{
		_enableStdReaction = reaction;
	}

	bool
	getEnableReSolution() const noexcept
	{
		return _enableReSolution;
	}

	virtual void
	setEnableReSolution(bool reso)
	{
		_enableReSolution = reso;
	}

	bool
	getEnableNucleation() const noexcept
	{
		return _enableNucleation;
	}

	virtual void
	setEnableNucleation(bool nuc)
	{
		_enableNucleation = nuc;
	}

	bool
	getEnableSink() const noexcept
	{
		return _enableSink;
	}

	virtual void
	setEnableSink(bool sink)
	{
		_enableSink = sink;
	}

	bool
	getEnableTrapMutation() const noexcept
	{
		return _enableTrapMutation;
	}

	virtual void
	setEnableTrapMutation(bool tm)
	{
		_enableTrapMutation = tm;
	}

	bool
	getEnableAttenuation() const noexcept
	{
		return _enableAttenuation;
	}

	virtual void
	setEnableAttenuation(bool enable)
	{
		_enableAttenuation = enable;
	}

	bool
	getEnableConstantReaction() const noexcept
	{
		return _enableConstantReaction;
	}

	virtual void
	setEnableConstantReaction(bool enable)
	{
		_enableConstantReaction = enable;
	}

	bool
	getEnableReducedJacobian() const noexcept
	{
		return _enableReducedJacobian;
	}

	virtual void
	setEnableReducedJacobian(bool reduced)
	{
		_enableReducedJacobian = reduced;
	}

	bool
	getEnableReadRates() const noexcept
	{
		return _enableReadRates;
	}

	virtual void
	setEnableReadRates(bool read)
	{
		_enableReadRates = read;
	}

	IndexType
	getGridSize() const noexcept
	{
		return _gridSize;
	}

	virtual void
	setGridSize(IndexType gridSize) = 0;

	/**
	 * @brief Takes a vector of temperatures and associated depths along X and
	 * updates the diffusion coefficients and rates accordingly.
	 */
	virtual void
	setTemperatures(const std::vector<double>& gridTemperatures,
		const std::vector<double>& gridDepths) = 0;

	/**
	 * @brief To update time dependent rates.
	 */
	virtual void
	setTime(double time) = 0;

	/**
	 * @brief Copies tile and cluster data from device to host.
	 */
	virtual void
	syncClusterDataOnHost() = 0;

	virtual IndexType
	findClusterId(const std::vector<AmountType>& composition) = 0;

	virtual ClusterCommon<plsm::HostMemSpace>
	getClusterCommon(IndexType clusterId) = 0;

	virtual ClusterCommon<plsm::HostMemSpace>
	getSingleVacancy() = 0;

	virtual IndexType
	getLargestClusterId() = 0;

	/**
	 * @brief Returns an object representing the the bounds of each
	 * cluster in each dimension of the phase space.
	 */
	virtual Bounds
	getAllClusterBounds() = 0;

	/**
	 * @brief Returns an object representing the the bounds of each
	 * cluster in each dimension of the phase space.
	 */
	virtual MomentIdMap
	getAllMomentIdInfo() = 0;

	/**
	 * @brief Return a string of cluster name in ID order.
	 */
	virtual std::string
	getHeaderString() = 0;

	/**
	 * @brief Computes the map between the different cluster bounds and moment
	 * IDs.
	 */
	virtual void initializeClusterMap(
		BoundVector, MomentIdMapVector, MomentIdMap) = 0;

	/**
	 * @brief Initialize reactions in the case it was not already
	 * done in the constructor.
	 */
	virtual void
	initializeReactions() = 0;

	/**
	 * @brief Set the rates for constant reactions
	 */
	virtual void setConstantRates(RateVector) = 0;

	/**
	 * @brief Set the connectivities for constant reactions
	 */
	virtual void setConstantConnectivities(ConnectivitiesVector) = 0;

	virtual PhaseSpace
	getPhaseSpace() = 0;

	/**
	 * @brief Updates the fluxes view with the rates from all the
	 * reactions at this grid point, the fluxes are used by the RHS function.
	 */
	virtual void
	computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex = 0, double surfaceDepth = 0.0,
		double spacing = 0.0) = 0;

	/**
	 * @brief Updates the values view with the rates from all the
	 * reactions at this grid point, they are used by the RHS Jacobian.
	 */
	virtual void
	computeAllPartials(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex = 0,
		double surfaceDepth = 0.0, double spacing = 0.0) = 0;

	/**
	 * @brief Updates the rates view with the rates from all the
	 * reactions at this grid point, this is for multiple instances use.
	 */
	virtual void
	computeConstantRates(ConcentrationsView concentrations, RatesView rates,
		IndexType subId, IndexType gridIndex = 0, double surfaceDepth = 0.0,
		double spacing = 0.0) = 0;

	/**
	 * @brief Updates the rates view with the rates from all the
	 * reactions at this grid point, this is for multiple instances use.
	 */
	virtual void
	getConstantConnectivities(ConnectivitiesView conns, IndexType subId) = 0;

	/**
	 * @brief Returns the largest computed rate.
	 */
	virtual double
	getLargestRate() = 0;

	/**
	 * @brief Returns the sum of rates at this grid point from reactions
	 * where the given cluster Id is on the left side of the reaction (either
	 * a reactant or dissociating).
	 */
	virtual double
	getLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex) = 0;

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 * Also populates the inverse map.
	 *
	 * @param fillMap Connectivity map.
	 * @return The total number of partials.
	 */
	virtual IndexType
	getDiagonalFill(SparseFillMap& fillMap) = 0;

	virtual double
	getTotalConcentration(ConcentrationsView concentrations, SpeciesId species,
		AmountType minSize = 0) = 0;

	virtual double
	getTotalRadiusConcentration(ConcentrationsView concentrations,
		SpeciesId species, AmountType minSize = 0) = 0;

	virtual double
	getTotalAtomConcentration(ConcentrationsView concentrations,
		SpeciesId species, AmountType minSize = 0) = 0;

	/**
	 * @brief Computes the diffusion flux exiting from this grid point.
	 *
	 * @param gridPointSolution The array of local solution, indexed with
	 * cluster Ids
	 * @param factor The geometric factor
	 * @param diffusingIds The vector of cluster Ids that we want to compute the
	 * diffusion flux for
	 * @param fluxes The vector of fluxes to update
	 * @param gridIndex The grid point location
	 */
	virtual void
	updateOutgoingDiffFluxes(double* gridPointSolution, double factor,
		std::vector<IndexType> diffusingIds, std::vector<double>& fluxes,
		IndexType gridIndex) = 0;

	/**
	 * @brief Computes the advection flux exiting from this grid point.
	 *
	 * @param gridPointSolution The array of local solution, indexed with
	 * cluster Ids
	 * @param factor The geometric factor
	 * @param advectingIds The vector of cluster Ids that we want to compute the
	 * advection flux for
	 * @param sinkStrengths The vector of sink strengths corresponding to each
	 * cluster advection
	 * @param fluxes The vector of fluxes to update
	 * @param gridIndex The grid point location
	 */
	virtual void
	updateOutgoingAdvecFluxes(double* gridPointSolution, double factor,
		std::vector<IndexType> advectingIds, std::vector<double> sinkStrengths,
		std::vector<double>& fluxes, IndexType gridIndex) = 0;

protected:
	std::string _material;
	double _latticeParameter{};
	double _atomicVolume{};
	double _interstitialBias{};
	double _impurityRadius{};
	double _fissionRate{};
	bool _enableStdReaction{};
	bool _enableReSolution{};
	bool _enableNucleation{};
	bool _enableSink{};
	bool _enableTrapMutation{};
	bool _enableAttenuation{};
	bool _enableConstantReaction{};
	bool _enableReducedJacobian{};
	bool _enableReadRates{};

	IndexType _gridSize{};
	IndexType _numDOFs{};
	IndexType _numClusters{};
};
} // namespace network
} // namespace core
} // namespace xolotl
