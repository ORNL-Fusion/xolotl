#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
struct DesorptionInitializer
{
	using AmountType = detail::CompositionAmountType;

	void
	set(AmountType dSize, double dPortion) noexcept
	{
		size = dSize;
		portion = dPortion;
	}

	AmountType size{};
	double portion{};
};

struct Desorption
{
	using AmountType = detail::CompositionAmountType;
	using IndexType = detail::ReactionNetworkIndexType;

	Desorption() = default;

	Desorption(DesorptionInitializer dInit, IndexType dId) :
		size{dInit.size},
		portion{dInit.portion},
		id{dId}
	{
	}

	AmountType size{};
	double portion{};
	IndexType id{detail::invalidNetworkIndex};
};

template <typename TClusterDataParent>
struct TrapMutationClusterData
{
	using AmountType = typename TClusterDataParent::AmountType;
	template <typename TData>
	using View = typename TClusterDataParent::template View<TData>;

	TrapMutationClusterData() = default;

	template <typename TClusterDataOther>
	KOKKOS_INLINE_FUNCTION
	TrapMutationClusterData(
		const TrapMutationClusterData<TClusterDataOther>& other) :
		desorption(other.desorption),
		currentDesorpLeftSideRate(other.currentDesorpLeftSideRate),
		currentDisappearingRate(other.currentDisappearingRate),
		tmDepths(other.tmDepths),
		tmVSizes(other.tmVSizes),
		tmEnabled(other.tmEnabled)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept;

	void
	initialize();

	View<Desorption> desorption;
	View<double> currentDesorpLeftSideRate;
	View<double> currentDisappearingRate;
	View<double[7]> tmDepths;
	View<AmountType[7]> tmVSizes;
	View<bool[7]> tmEnabled;
};

template <typename TP1, typename TP2>
inline void
deepCopy(TrapMutationClusterData<TP1> to, TrapMutationClusterData<TP2> from);
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
