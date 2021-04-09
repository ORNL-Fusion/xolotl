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
	View<double[7]> tmDepths; // should be DualView
	View<AmountType[7]> tmVSizes; // may only be needed at initialization
	View<bool[7]> tmEnabled;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
