#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
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
	}

	View<Desorption> desorption;
	View<double> currentDesorpLeftSideRate;
	View<double> currentDisappearingRate;
	View<double[7]> tmDepths; // should be DualView
	View<AmountType[7]> tmVSizes; // may only be needed at initialization
	View<bool[7]> tmEnabled;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
