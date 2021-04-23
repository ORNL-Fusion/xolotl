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
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = 0;

		ret += desorption.required_allocation_size();
		ret += currentDesorpLeftSideRate.required_allocation_size();
		ret += currentDisappearingRate.required_allocation_size();
		ret += tmDepths.required_allocation_size();
		ret += tmVSizes.required_allocation_size();
		ret += tmEnabled.required_allocation_size();

		return ret;
	}

	void
	initialize()
	{
		currentDesorpLeftSideRate =
			View<double>("Current Desorption Left Side Rate");

		currentDisappearingRate =
			View<double>("Current Trap Mutation Disappearing Rate");
		auto mirror = create_mirror_view(currentDisappearingRate);
		mirror() = 1.0;
		deep_copy(currentDisappearingRate, mirror);

		desorption = View<Desorption>("Desorption");
		tmDepths = View<double[7]>("Trap-mutation depths");
		tmVSizes = View<AmountType[7]>("Trap-mutation vacancy sizes");

		tmEnabled = View<bool[7]>("Trap-mutation enabled helium sizes");
	}

	View<Desorption> desorption;
	View<double> currentDesorpLeftSideRate;
	View<double> currentDisappearingRate;
	View<double[7]> tmDepths;
	View<AmountType[7]> tmVSizes;
	View<bool[7]> tmEnabled;
};

template <typename TP1, typename TP2>
inline void
deepCopy(TrapMutationClusterData<TP1> to, TrapMutationClusterData<TP2> from)
{
	if (!from.desorption.is_allocated()) {
		return;
	}

	if (!to.desorption.is_allocated()) {
		to.desorption = create_mirror_view(from.desorption);
		to.currentDesorpLeftSideRate =
			create_mirror_view(from.currentDesorpLeftSideRate);
		to.currentDisappearingRate =
			create_mirror_view(from.currentDisappearingRate);
		to.tmDepths = create_mirror_view(from.tmDepths);
		to.tmVSizes = create_mirror_view(from.tmVSizes);
		to.tmEnabled = create_mirror_view(from.tmEnabled);
	}

	deep_copy(to.desorption, from.desorption);
	deep_copy(to.currentDesorpLeftSideRate, from.currentDesorpLeftSideRate);
	deep_copy(to.currentDisappearingRate, from.currentDisappearingRate);
	deep_copy(to.tmDepths, from.tmDepths);
	deep_copy(to.tmVSizes, from.tmVSizes);
	deep_copy(to.tmEnabled, from.tmEnabled);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
