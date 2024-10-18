#pragma once

#include <string>

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename MemSpace>
class ClusterCommon;

template <typename TNetwork, typename MemSpace>
class Cluster;

namespace detail
{
template <typename TData, typename MemSpace>
struct ViewTypeHelper
{
	using DeviceView = Kokkos::View<TData, plsm::DeviceMemSpace>;
	using HostView = Kokkos::View<TData, plsm::HostMemSpace>;
	using ViewType = std::conditional_t<
		std::is_same_v<plsm::HostMemSpace, plsm::DeviceMemSpace>, HostView,
		std::conditional_t<std::is_same_v<MemSpace, plsm::DeviceMemSpace>,
			DeviceView, typename DeviceView::HostMirror>>;
};

template <typename TData, typename MemSpace>
using ViewType = typename ViewTypeHelper<TData, MemSpace>::ViewType;

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

template <typename TNetwork, typename MemSpace>
struct ClusterDataExtra
{
	using IndexType = detail::ReactionNetworkIndexType;
	static_assert(Kokkos::is_memory_space<MemSpace>{});

	ClusterDataExtra() = default;

	template <typename MS>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<TNetwork, MS>&)
	{
	}

	template <typename MS>
	void
	deepCopy([[maybe_unused]] const ClusterDataExtra<TNetwork, MS>& data)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return 0;
	}

	void
	setGridSize(IndexType numClusters, IndexType gridSize) {};
};

/**
 * @brief Structure for physical properties and clusters,
 * independent of the network type.
 *
 * @tparam MemSpace plsm::HostMemSpace or plsm::DeviceMemSpace
 */
template <typename MemSpace>
struct ClusterDataCommon
{
	static_assert(Kokkos::is_memory_space<MemSpace>{});

	template <typename>
	friend class ClusterDataCommon;

	template <typename TData>
	using View = ViewType<TData, MemSpace>;

	using ClusterType = ClusterCommon<MemSpace>;
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
		_intVals(data._intVals),
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
		AV_HELIUM_CONTENT,
		AV_VACANCY_CONTENT,
		AV_RADIUS,
		DEPTH,
		TAU_BURSTING,
		F_BURSTING,
		TEMPERATURE,
		I_FORMATION,
		XE_FORMATION,
		V_FORMATION,
		V2_FORMATION,
		NUM_FLOAT_VALS
	};

	enum IntValsIndex : int
	{
		BUBBLE = 0,
		AV_HELIUM,
		AV_VACANCY,
		TRANSITION_SIZE,
		NUM_INT_VALS
	};

	enum BoolValsIndex : int
	{
		STD_REACTION = 0,
		RESOLUTION,
		NUCLEATION,
		SINK,
		TRAP_MUTATION,
		BURST,
		LARGE_BUBBLE,
		READ_RATES,
		CONSTANT_REACTION,
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
	double
	bubbleAvHe() const
	{
		return _floatVals[AV_HELIUM_CONTENT];
	}

	void
	setBubbleAvHe(double val)
	{
		setVal(_floatVals, AV_HELIUM_CONTENT, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getIFormationEnergy() const
	{
		return _floatVals[I_FORMATION];
	}

	void
	setIFormationEnergy(double val)
	{
		setVal(_floatVals, I_FORMATION, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	bubbleAvV() const
	{
		return _floatVals[AV_VACANCY_CONTENT];
	}

	void
	setBubbleAvV(double val)
	{
		setVal(_floatVals, AV_VACANCY_CONTENT, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getVFormationEnergy() const
	{
		return _floatVals[V_FORMATION];
	}

	void
	setVFormationEnergy(double val)
	{
		setVal(_floatVals, V_FORMATION, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	bubbleAvRadius() const
	{
		return _floatVals[AV_RADIUS];
	}

	void
	setBubbleAvRadius(double val)
	{
		setVal(_floatVals, AV_RADIUS, val);
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	bubbleId() const
	{
		return _intVals[BUBBLE];
	}

	void
	setBubbleId(IndexType val)
	{
		setVal(_intVals, BUBBLE, val);
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	bubbleAvHeId() const
	{
		return _intVals[AV_HELIUM];
	}

	void
	setBubbleAvHeId(IndexType val)
	{
		setVal(_intVals, AV_HELIUM, val);
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	bubbleAvVId() const
	{
		return _intVals[AV_VACANCY];
	}

	void
	setBubbleAvVId(IndexType val)
	{
		setVal(_intVals, AV_VACANCY, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getV2FormationEnergy() const
	{
		return _floatVals[V2_FORMATION];
	}

	void
	setV2FormationEnergy(double val)
	{
		setVal(_floatVals, V2_FORMATION, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getDepth() const
	{
		return _floatVals[DEPTH];
	}

	void
	setDepth(double val)
	{
		setVal(_floatVals, DEPTH, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getTauBursting() const
	{
		return _floatVals[TAU_BURSTING];
	}

	void
	setTauBursting(double val)
	{
		setVal(_floatVals, TAU_BURSTING, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getFBursting() const
	{
		return _floatVals[F_BURSTING];
	}

	void
	setFBursting(double val)
	{
		setVal(_floatVals, F_BURSTING, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getTemperature() const
	{
		return _floatVals[TEMPERATURE];
	}

	void
	setTemperature(double val)
	{
		setVal(_floatVals, TEMPERATURE, val);
	}

	KOKKOS_INLINE_FUNCTION
	double
	getXeFormationEnergy() const
	{
		return _floatVals[XE_FORMATION];
	}

	void
	setXeFormationEnergy(double val)
	{
		setVal(_floatVals, XE_FORMATION, val);
	}

	int
	transitionSize() const
	{
		return _intVals[TRANSITION_SIZE];
	}

	void
	setTransitionSize(int val)
	{
		setVal(_intVals, TRANSITION_SIZE, val);
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

	KOKKOS_INLINE_FUNCTION
	bool
	enableBurst() const
	{
		return _boolVals[BURST];
	}

	void
	setEnableBurst(bool val)
	{
		setVal(_boolVals, BURST, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableLargeBubble() const
	{
		return _boolVals[LARGE_BUBBLE];
	}

	void
	setEnableLargeBubble(bool val)
	{
		setVal(_boolVals, LARGE_BUBBLE, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableReadRates() const
	{
		return _boolVals[READ_RATES];
	}

	void
	setEnableReadRates(bool val)
	{
		setVal(_boolVals, READ_RATES, val);
	}

	KOKKOS_INLINE_FUNCTION
	bool
	enableConstantReaction() const
	{
		return _boolVals[CONSTANT_REACTION];
	}

	void
	setEnableConstantReaction(bool val)
	{
		setVal(_boolVals, CONSTANT_REACTION, val);
	}

private:
	View<double[NUM_FLOAT_VALS]> _floatVals;
	View<IndexType[NUM_INT_VALS]> _intVals;
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
 * @tparam MemSpace plsm::HostMemSpace or plsm::DeviceMemSpace
 */
template <typename TNetwork, typename MemSpace>
struct ClusterData : ClusterDataCommon<MemSpace>
{
	static_assert(Kokkos::is_memory_space<MemSpace>{});

private:
	using Traits = ReactionNetworkTraits<TNetwork>;
	using Types = detail::ReactionNetworkTypes<TNetwork>;
	using Props = detail::ReactionNetworkProperties<TNetwork>;
	static constexpr auto nMomentIds = Props::numSpeciesNoI;

public:
	using Superclass = ClusterDataCommon<MemSpace>;
	using ClusterGenerator = typename Traits::ClusterGenerator;
	using ClusterUpdater = typename Types::ClusterUpdater;
	using Subpaving =
		plsm::MemSpaceSubpaving<MemSpace, typename Types::Subpaving>;
	using TilesView = typename Subpaving::TilesView;
	using ClusterType = Cluster<TNetwork, MemSpace>;
	using IndexType = typename Types::IndexType;

	template <typename TData>
	using View = typename Superclass::template View<TData>;

	ClusterData() = default;

	ClusterData(const TilesView& tiles_, IndexType numClusters_,
		IndexType gridSize_ = 0);

	explicit ClusterData(Subpaving& subpaving, IndexType gridSize_ = 0);

	template <typename TClusterData>
	KOKKOS_INLINE_FUNCTION
	ClusterData(const TClusterData& data) :
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

	void
	updateDiffusionCoefficients();

	IndexType
	defineMomentIds();

	TilesView tiles;
	View<IndexType* [nMomentIds]> momentIds;
	ClusterDataExtra<TNetwork, MemSpace> extraData;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
