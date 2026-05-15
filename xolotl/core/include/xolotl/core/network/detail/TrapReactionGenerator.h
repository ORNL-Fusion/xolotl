#pragma once

#include <xolotl/core/network/detail/ReactionGenerator.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief Inherits from ReactionGeneratorBase and additionally fills collections
 * trap reactions depending on the subpaving.
 *
 * @tparam TBase The templated base class.
 */
template <typename TBase>
class TrapReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using TrapReactionType = typename NetworkTraits::TrapReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	TrapReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addTrapReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addTrapReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<TrapReactionType*>
	getTrapReactions() const
	{
		return _trapReactions;
	}

	IndexType
	getNumberOfTrapReactions() const
	{
		return _trapReactions.size();
	}

private:
	IndexView _clusterTrapReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numTrapReactions{};

	IndexView _trapCrsRowMap;
	ClusterSetSubView _trapCrsClusterSets;

	Kokkos::View<TrapReactionType*> _trapReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<
		std::is_base_of_v<TrapReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = TrapReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
