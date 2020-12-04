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
 * sink reactions depending on the subpaving.
 *
 * @tparam TBase The templated base class.
 */
template <typename TBase>
class SinkReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using SinkReactionType = typename NetworkTraits::SinkReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	SinkReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addSinkReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addSinkReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<SinkReactionType*>
	getSinkReactions() const
	{
		return _sinkReactions;
	}

	IndexType
	getNumberOfSinkReactions() const
	{
		return _sinkReactions.size();
	}

private:
	IndexView _clusterSinkReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numSinkReactions{};

	IndexView _sinkCrsRowMap;
	ClusterSetSubView _sinkCrsClusterSets;

	Kokkos::View<SinkReactionType*> _sinkReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<
		std::is_base_of<SinkReaction<TNetwork, TReaction>, TReaction>::value>>
{
	using Type = SinkReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
