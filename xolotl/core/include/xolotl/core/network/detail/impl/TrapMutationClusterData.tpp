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
inline std::uint64_t
TrapMutationClusterData<TClusterDataParent>::getDeviceMemorySize()
	const noexcept
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

template <typename TClusterDataParent>
inline void
TrapMutationClusterData<TClusterDataParent>::initialize()
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
