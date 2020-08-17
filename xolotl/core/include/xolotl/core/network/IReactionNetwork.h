#pragma once

#include <unordered_map>

#include <Kokkos_Core.hpp>
#include <Kokkos_Crs.hpp>

#include <xolotl/core/network/Cluster.h>
#include <xolotl/core/network/ClusterData.h>
#include <xolotl/core/network/SpeciesId.h>
#include <xolotl/core/network/detail/ReactionData.h>

namespace xolotl
{
namespace core
{
namespace network
{
class IReactionNetwork
{
public:
	using IndexType = detail::ReactionNetworkIndexType;
	using AmountType = detail::CompositionAmountType;
	using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
	using OwnedConcentrationsView = Kokkos::View<double*>;
	using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
	using OwnedFluxesView = Kokkos::View<double*>;
	using Connectivity = detail::ClusterConnectivity<>;
	using SparseFillMap = std::unordered_map<int, std::vector<int>>;
	using Bounds = std::vector<std::vector<AmountType>>;
	using PhaseSpace = std::vector<const char*>;

	static constexpr IndexType
	invalidIndex() noexcept
	{
		return detail::InvalidIndex::value;
	}

	static constexpr AmountType
	invalidAmount() noexcept
	{
		return detail::InvalidSpeciesAmount::value;
	}

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
	setFissionRate(double rate) = 0;

	virtual void
	setZeta(double z) = 0;

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableStdReaction() const noexcept
	{
		return _enableStdReaction;
	}

	virtual void
	setEnableStdReaction(bool reaction) = 0;

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableReSolution() const noexcept
	{
		return _enableReSolution;
	}

	virtual void
	setEnableReSolution(bool reso) = 0;

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableNucleation() const noexcept
	{
		return _enableNucleation;
	}

	virtual void
	setEnableNucleation(bool reso) = 0;

	KOKKOS_INLINE_FUNCTION
	IndexType
	getGridSize() const noexcept
	{
		return _gridSize;
	}

	virtual void
	setGridSize(IndexType gridSize) = 0;

	virtual void
	setTemperatures(const std::vector<double>& gridTemperatures) = 0;

	virtual void
	syncClusterDataOnHost() = 0;

    virtual IndexType
    findClusterId(const std::vector<AmountType>& composition) = 0;

	virtual ClusterCommon<plsm::OnHost>
	getClusterCommon(IndexType clusterId) const = 0;

	virtual ClusterCommon<plsm::OnHost>
	getSingleVacancy() = 0;

	virtual IndexType
	getLargestClusterId() = 0;

	virtual Bounds
	getAllClusterBounds() = 0;

	virtual PhaseSpace
	getPhaseSpace() = 0;

	virtual void
	computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex) = 0;

	virtual void
	computeAllPartials(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex) = 0;

	virtual double
	getLargestRate() = 0;

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
	getTotalAtomConcentration(ConcentrationsView concentrations,
		SpeciesId species, AmountType minSize = 0) = 0;

protected:
	double _latticeParameter{};
	double _atomicVolume{};
	double _interstitialBias{};
	double _impurityRadius{};
	double _fissionRate{};
	bool _enableStdReaction{};
	bool _enableReSolution{};
	bool _enableNucleation{};

	IndexType _gridSize{};
	IndexType _numDOFs{};
	IndexType _numClusters{};
};
} // namespace network
} // namespace core
} // namespace xolotl
