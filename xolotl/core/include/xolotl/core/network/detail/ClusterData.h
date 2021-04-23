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
	template <typename TData>
	using View = ViewConvert<ViewType<TData, PlsmContext>>;

	using ClusterType = ClusterCommon<PlsmContext>;
	using IndexType = detail::ReactionNetworkIndexType;
	using AmountType = detail::CompositionAmountType;

	ClusterDataCommon() = default;

	explicit ClusterDataCommon(IndexType numClusters_, IndexType gridSize_ = 0);

	template <typename TClusterDataCommon>
	KOKKOS_INLINE_FUNCTION
	ClusterDataCommon(const TClusterDataCommon& data) :
		numClusters(data.numClusters),
		gridSize(data.gridSize),
		_floatVals(data._floatVals),
		atomicVolume(data.atomicVolume),
		latticeParameter(data.latticeParameter),
		fissionRate(data.fissionRate),
		zeta(data.zeta),
		_boolVals(data._boolVals),
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
	getDeviceMemorySize() const noexcept;

	ClusterType
	getCluster(IndexType clusterId) const noexcept
	{
		return ClusterType(*this, clusterId);
	}

	void
	setGridSize(IndexType gridSize_);

	IndexType numClusters{};
	IndexType gridSize{};
	View<double[4]> _floatVals;
	Unmanaged<View<double>> atomicVolume;
	Unmanaged<View<double>> latticeParameter;
	Unmanaged<View<double>> fissionRate;
	Unmanaged<View<double>> zeta;
	View<bool[5]> _boolVals;
	Unmanaged<View<bool>> enableStdReaction;
	Unmanaged<View<bool>> enableReSolution;
	Unmanaged<View<bool>> enableNucleation;
	Unmanaged<View<bool>> enableSink;
	Unmanaged<View<bool>> enableTrapMutation;

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
	using Traits = ReactionNetworkTraits<TNetwork>;
	using Types = detail::ReactionNetworkTypes<TNetwork>;
	using Props = detail::ReactionNetworkProperties<TNetwork>;
	static constexpr auto nMomentIds = Props::numSpeciesNoI;

public:
	using Superclass = ClusterDataCommon<PlsmContext, ViewConvert>;
	using ClusterGenerator = typename Traits::ClusterGenerator;
	using Subpaving = typename Types::Subpaving;
	using TilesView =
		ViewConvert<typename Subpaving::template TilesView<PlsmContext>>;
	using ClusterType = Cluster<TNetwork, PlsmContext>;
	using IndexType = typename Types::IndexType;

	template <typename TData>
	using View = typename Superclass::template View<TData>;

	ClusterDataImpl() = default;

	ClusterDataImpl(const TilesView& tiles_, IndexType numClusters_,
		IndexType gridSize_ = 0);

	explicit ClusterDataImpl(Subpaving& subpaving, IndexType gridSize_ = 0);

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
	getDeviceMemorySize() const noexcept;

	KOKKOS_INLINE_FUNCTION
	ClusterType
	getCluster(IndexType clusterId) const noexcept
	{
		return ClusterType(*this, clusterId);
	}

	void
	generate(const ClusterGenerator& generator, double latticeParameter,
		double interstitialBias, double impurityRadius);

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
deepCopy(
	ClusterDataImpl<TN1, PC1, VC1> to, ClusterDataImpl<TN2, PC2, VC2> from);
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
