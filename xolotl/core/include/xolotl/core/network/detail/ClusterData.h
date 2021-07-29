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

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
struct ClusterDataExtra
{
	ClusterDataExtra() = default;

	template <typename PC, template <typename> typename VC>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<TNetwork, PC, VC>&)
	{
	}

	template <typename PC, template <typename> typename VC>
	void
	deepCopy([[maybe_unused]] const ClusterDataExtra<TNetwork, PC, VC>& data)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return 0;
	}
};

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
	template <typename, template <typename> typename>
	friend class ClusterDataCommon;

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
		_boolVals(data._boolVals),
		temperature(data.temperature),
		reactionRadius(data.reactionRadius),
		formationEnergy(data.formationEnergy),
		migrationEnergy(data.migrationEnergy),
		diffusionFactor(data.diffusionFactor),
		diffusionCoefficient(data.diffusionCoefficient)
	{
	}

	template <typename TClusterDataCommon>
	void
	deepCopy(const TClusterDataCommon& data);

	std::uint64_t
	getDeviceMemorySize() const noexcept;

	ClusterType
	getCluster(IndexType clusterId) const noexcept
	{
		return ClusterType(this, clusterId);
	}

	ClusterType
	getClusterCommon(IndexType clusterId) const noexcept
	{
		return getCluster(clusterId);
	}

	void
	setGridSize(IndexType gridSize_);

private:
	enum FloatValsIndex : int
	{
		ATOMIC_VOLUME = 0,
		LATTICE_PARAM,
		FISSION_RATE,
		ZETA,
		NUM_FLOAT_VALS
	};

	enum BoolValsIndex : int
	{
		STD_REACTION = 0,
		RESOLUTION,
		NUCLEATION,
		SINK,
		TRAP_MUTATION,
		NUM_BOOL_VALS
	};

	template <typename TView, typename TVal>
	void
	setVal(TView view, int index, TVal value)
	{
		auto sub = subview(view, index);
		auto mir = create_mirror_view(sub);
		mir() = value;
		deep_copy(sub, mir);
	}

public:
	KOKKOS_INLINE_FUNCTION
	double
	atomicVolume() const
	{
		return _floatVals[ATOMIC_VOLUME];
	}

	void
	setAtomicVolume(double val)
	{
		setVal(_floatVals, ATOMIC_VOLUME, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	latticeParameter() const
	{
		return _floatVals[LATTICE_PARAM];
	}

	void
	setLatticeParameter(double val)
	{
		setVal(_floatVals, LATTICE_PARAM, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	fissionRate() const
	{
		return _floatVals[FISSION_RATE];
	}

	void
	setFissionRate(double val)
	{
		setVal(_floatVals, FISSION_RATE, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	zeta() const
	{
		return _floatVals[ZETA];
	}

	void
	setZeta(double val)
	{
		setVal(_floatVals, ZETA, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableStdReaction() const
	{
		return _boolVals[STD_REACTION];
	}

	void
	setEnableStdReaction(bool val)
	{
		setVal(_boolVals, STD_REACTION, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableReSolution() const
	{
		return _boolVals[RESOLUTION];
	}

	void
	setEnableReSolution(bool val)
	{
		setVal(_boolVals, RESOLUTION, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableNucleation() const
	{
		return _boolVals[NUCLEATION];
	}

	void
	setEnableNucleation(bool val)
	{
		setVal(_boolVals, NUCLEATION, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableSink() const
	{
		return _boolVals[SINK];
	}

	void
	setEnableSink(bool val)
	{
		setVal(_boolVals, SINK, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableTrapMutation() const
	{
		return _boolVals[TRAP_MUTATION];
	}

	void
	setEnableTrapMutation(bool val)
	{
		setVal(_boolVals, TRAP_MUTATION, val);
	}

private:
	View<double[NUM_FLOAT_VALS]> _floatVals;
	View<bool[NUM_BOOL_VALS]> _boolVals;

public:
	IndexType numClusters{};
	IndexType gridSize{};

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

	template <typename TClusterData>
	void
	deepCopy(const TClusterData& data);

	std::uint64_t
	getDeviceMemorySize() const noexcept;

	KOKKOS_INLINE_FUNCTION
	ClusterType
	getCluster(IndexType clusterId) const noexcept
	{
		return ClusterType(this, clusterId);
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
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
