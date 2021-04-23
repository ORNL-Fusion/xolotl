#pragma once

#include <string>

#include <Kokkos_View.hpp>

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/MemorySpace.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename PlsmContext>
class ClusterCommon;

template <typename TNetwork, typename PlsmContext>
class Cluster;

namespace detail
{
template <typename PlsmContext>
struct ContextLabelHelper;

template <>
struct ContextLabelHelper<plsm::OnHost>
{
	static constexpr char label = 'H';
};

template <>
struct ContextLabelHelper<plsm::OnDevice>
{
	static constexpr char label = 'D';
};

template <typename PlsmContext>
constexpr char contextLabel = ContextLabelHelper<PlsmContext>::label;

inline std::string
labelStr(const char labelChar)
{
	return std::string("_") + labelChar;
}

template <typename TData, typename PlsmContext>
struct ViewTypeHelper;

template <typename TData>
struct ViewTypeHelper<TData, plsm::OnHost>
{
	using DeviceView = Kokkos::View<TData, DefaultMemorySpace>;
	using ViewType = typename DeviceView::HostMirror;
};

template <typename TData>
struct ViewTypeHelper<TData, plsm::OnDevice>
{
	using ViewType = Kokkos::View<TData, DefaultMemorySpace>;
};

template <typename TData, typename PlsmContext>
using ViewType = typename ViewTypeHelper<TData, PlsmContext>::ViewType;

template <typename TView>
struct UnmanagedHelper
{
	using Traits = typename TView::traits;
	using Type =
		Kokkos::View<typename Traits::data_type, typename Traits::array_layout,
			typename Traits::device_type, Kokkos::MemoryUnmanaged>;
};

template <typename TView>
using Unmanaged = typename UnmanagedHelper<TView>::Type;

template <typename TView>
using PassThru = TView;

/**
 * @brief Structure for physical properties and clusters,
 * independent of the network type.
 *
 * @tparam PlsmContext Host or Device
 */
template <typename PlsmContext,
	template <typename> typename ViewConvert = PassThru>
struct ClusterDataCommon
{
protected:
	static constexpr char label = contextLabel<PlsmContext>;

public:
	template <typename TData>
	using View = ViewConvert<ViewType<TData, PlsmContext>>;

	using ClusterType = ClusterCommon<PlsmContext>;
	using IndexType = detail::ReactionNetworkIndexType;
	using AmountType = detail::CompositionAmountType;

	ClusterDataCommon() = default;

	explicit ClusterDataCommon(
		IndexType numClusters_, IndexType gridSize_ = 0) :
		numClusters(numClusters_),
		gridSize(gridSize_),
		atomicVolume("Atomic Volume" + labelStr(label)),
		latticeParameter("Lattice Parameter" + labelStr(label)),
		fissionRate("Fission Rate" + labelStr(label)),
		zeta("Zeta" + labelStr(label)),
		enableStdReaction("Enable Std Reaction" + labelStr(label)),
		enableReSolution("Enable Re-Solution Process" + labelStr(label)),
		enableNucleation("Enable Nucleation Process" + labelStr(label)),
		enableSink("Enable Sink Process" + labelStr(label)),
		enableTrapMutation("Enable Trap Mutation Process" + labelStr(label)),
		temperature("Temperature" + labelStr(label), gridSize),
		reactionRadius("Reaction Radius" + labelStr(label), numClusters),
		formationEnergy("Formation Energy" + labelStr(label), numClusters),
		migrationEnergy("Migration Energy" + labelStr(label), numClusters),
		diffusionFactor("Diffusion Factor" + labelStr(label), numClusters),
		diffusionCoefficient(
			"Diffusion Coefficient" + labelStr(label), numClusters, gridSize)
	{
	}

	template <typename TClusterDataCommon>
	KOKKOS_INLINE_FUNCTION
	ClusterDataCommon(const TClusterDataCommon& data) :
		numClusters(data.numClusters),
		gridSize(data.gridSize),
		atomicVolume(data.atomicVolume),
		latticeParameter(data.latticeParameter),
		fissionRate(data.fissionRate),
		zeta(data.zeta),
		enableStdReaction(data.enableStdReaction),
		enableReSolution(data.enableReSolution),
		enableNucleation(data.enableNucleation),
		enableSink(data.enableSink),
		enableTrapMutation(data.enableTrapMutation),
		temperature(data.temperature),
		reactionRadius(data.reactionRadius),
		formationEnergy(data.formationEnergy),
		migrationEnergy(data.migrationEnergy),
		diffusionFactor(data.diffusionFactor),
		diffusionCoefficient(data.diffusionCoefficient)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = 0;

		ret += sizeof(numClusters);
		ret += sizeof(gridSize);
		ret += atomicVolume.required_allocation_size();
		ret += latticeParameter.required_allocation_size();
		ret += fissionRate.required_allocation_size();
		ret += zeta.required_allocation_size();
		ret += enableStdReaction.required_allocation_size();
		ret += enableReSolution.required_allocation_size();
		ret += enableNucleation.required_allocation_size();
		ret += enableSink.required_allocation_size();
		ret += enableTrapMutation.required_allocation_size();

		ret += temperature.required_allocation_size(temperature.size());
		ret += reactionRadius.required_allocation_size(reactionRadius.size());
		ret += formationEnergy.required_allocation_size(formationEnergy.size());
		ret += migrationEnergy.required_allocation_size(migrationEnergy.size());
		ret += diffusionFactor.required_allocation_size(diffusionFactor.size());
		ret += diffusionCoefficient.required_allocation_size(
			diffusionCoefficient.extent(0), diffusionCoefficient.extent(1));

		return ret;
	}

	void
	syncClusterDataOnHost()
	{
	}

	ClusterType
	getCluster(IndexType clusterId) const noexcept
	{
		return ClusterType(*this, clusterId);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getAtomicVolume() const
	{
		return atomicVolume(0);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getLatticeParameter() const
	{
		return latticeParameter(0);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getFissionRate() const
	{
		return fissionRate(0);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getZeta() const
	{
		return zeta(0);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableStdReaction() const
	{
		return enableStdReaction(0);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableReSolution() const
	{
		return enableReSolution(0);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableNucleation() const
	{
		return enableNucleation(0);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableSink() const
	{
		return enableSink(0);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	getEnableTrapMutation() const
	{
		return enableTrapMutation();
	}

	void
	setGridSize(IndexType gridSize_)
	{
		gridSize = gridSize_;
		temperature = View<double*>("Temperature" + labelStr(label), gridSize);
		diffusionCoefficient = View<double**>(
			"Diffusion Coefficient" + labelStr(label), numClusters, gridSize);
	}

	IndexType numClusters{};
	IndexType gridSize{};
	View<double[1]> atomicVolume;
	View<double[1]> latticeParameter;
	View<double[1]> fissionRate;
	View<double[1]> zeta;
	View<bool[1]> enableStdReaction;
	View<bool[1]> enableReSolution;
	View<bool[1]> enableNucleation;
	View<bool[1]> enableSink;
	View<bool> enableTrapMutation;

	View<double*> temperature;
	View<double*> reactionRadius;
	View<double*> formationEnergy;
	View<double*> migrationEnergy;
	View<double*> diffusionFactor;
	View<double**> diffusionCoefficient;
};

/**
 * @brief Structure for additional clusters properties that are
 * dependent on the network type (tiles and moments).
 *
 * @tparam TNetwork The network type
 * @tparam PlsmContext Host or Device
 */
template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
struct ClusterDataImpl : ClusterDataCommon<PlsmContext, ViewConvert>
{
private:
	using Types = detail::ReactionNetworkTypes<TNetwork>;
	using Props = detail::ReactionNetworkProperties<TNetwork>;
	static constexpr auto nMomentIds = Props::numSpeciesNoI;

public:
	using Superclass = ClusterDataCommon<PlsmContext, ViewConvert>;
	using Subpaving = typename Types::Subpaving;
	using TilesView =
		ViewConvert<typename Subpaving::template TilesView<PlsmContext>>;
	using ClusterType = Cluster<TNetwork, PlsmContext>;
	using IndexType = typename Types::IndexType;

	template <typename TData>
	using View = typename Superclass::template View<TData>;

	ClusterDataImpl() = default;

	ClusterDataImpl(const TilesView& tiles_, IndexType numClusters_,
		IndexType gridSize_ = 0) :
		Superclass(numClusters_, gridSize_),
		tiles(tiles_),
		momentIds("Moment Ids" + labelStr(this->label), numClusters_)
	{
	}

	explicit ClusterDataImpl(Subpaving& subpaving, IndexType gridSize_ = 0) :
		ClusterDataImpl(subpaving.getTiles(PlsmContext{}),
			subpaving.getNumberOfTiles(PlsmContext{}), gridSize_)
	{
	}

	template <typename TClusterData>
	KOKKOS_INLINE_FUNCTION
	ClusterDataImpl(const TClusterData& data) :
		Superclass(data),
		tiles(data.tiles),
		momentIds(data.momentIds),
		extraData(data.extraData)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = Superclass::getDeviceMemorySize();

		ret += tiles.required_allocation_size(tiles.size());
		ret += momentIds.required_allocation_size(momentIds.size());
		ret += extraData.getDeviceMemorySize();

		return ret;
	}

	KOKKOS_INLINE_FUNCTION
	ClusterType
	getCluster(IndexType clusterId) const noexcept
	{
		return ClusterType(*this, clusterId);
	}

	TilesView tiles;
	View<IndexType* [nMomentIds]> momentIds;
	ClusterDataExtra<TNetwork, PlsmContext, ViewConvert> extraData;
};

template <typename TNetwork, typename PlsmContext = plsm::OnDevice>
struct ClusterDataHelper
{
	using Type = ClusterDataImpl<TNetwork, PlsmContext, PassThru>;
};

template <typename PlsmContext>
using ClusterDataCommonRef = ClusterDataCommon<PlsmContext, Unmanaged>;

template <typename TNetwork, typename PlsmContext = plsm::OnDevice>
struct ClusterDataRefHelper
{
	using Type = ClusterDataImpl<TNetwork, PlsmContext, Unmanaged>;
};

template <typename TN1, typename PC1, template <typename> typename VC1,
	typename TN2, typename PC2, template <typename> typename VC2>
inline void
deepCopy(ClusterDataImpl<TN1, PC1, VC1> to, ClusterDataImpl<TN2, PC2, VC2> from)
{
	deep_copy(to.atomicVolume, from.atomicVolume);
	deep_copy(to.latticeParameter, from.latticeParameter);
	deep_copy(to.fissionRate, from.fissionRate);
	deep_copy(to.zeta, from.zeta);
	deep_copy(to.enableStdReaction, from.enableStdReaction);
	deep_copy(to.enableReSolution, from.enableReSolution);
	deep_copy(to.enableNucleation, from.enableNucleation);
	deep_copy(to.enableSink, from.enableSink);
	deep_copy(to.enableTrapMutation, from.enableTrapMutation);
	deep_copy(to.temperature, from.temperature);
	deep_copy(to.reactionRadius, from.reactionRadius);
	deep_copy(to.formationEnergy, from.formationEnergy);
	deep_copy(to.migrationEnergy, from.migrationEnergy);
	deep_copy(to.diffusionFactor, from.diffusionFactor);
	deep_copy(to.diffusionCoefficient, from.diffusionCoefficient);

	// NOTE: Intentionally omitting tiles assuming that was part of construction
	deep_copy(to.momentIds, from.momentIds);

	deepCopy(to.extraData, from.extraData);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
