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
template <typename TClusterDataOther>
inline void
TrapMutationClusterData<TClusterDataParent>::deepCopy(
	const TrapMutationClusterData<TClusterDataOther>& from)
{
	if (!from.desorption.is_allocated()) {
		return;
	}

	if (!desorption.is_allocated()) {
		desorption = create_mirror_view(from.desorption);
		currentDesorpLeftSideRate =
			create_mirror_view(from.currentDesorpLeftSideRate);
		currentDisappearingRate =
			create_mirror_view(from.currentDisappearingRate);
		tmDepths = create_mirror_view(from.tmDepths);
		tmVSizes = create_mirror_view(from.tmVSizes);
		tmEnabled = create_mirror_view(from.tmEnabled);
	}

	deep_copy(desorption, from.desorption);
	deep_copy(currentDesorpLeftSideRate, from.currentDesorpLeftSideRate);
	deep_copy(currentDisappearingRate, from.currentDisappearingRate);
	deep_copy(tmDepths, from.tmDepths);
	deep_copy(tmVSizes, from.tmVSizes);
	deep_copy(tmEnabled, from.tmEnabled);
}

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
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
