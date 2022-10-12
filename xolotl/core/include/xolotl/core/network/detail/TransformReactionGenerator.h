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
 * transform reactions depending on the subpaving.
 *
 * @tparam TBase The templated base class.
 */
template <typename TBase>
class TransformReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using TransformReactionType = typename NetworkTraits::TransformReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	TransformReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addTransformReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addTransformReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<TransformReactionType*>
	getTransformReactions() const
	{
		return _transformReactions;
	}

	IndexType
	getNumberOfTransformReactions() const
	{
		return _transformReactions.size();
	}

private:
	IndexView _clusterTransformReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numTransformReactions{};

	IndexView _transformCrsRowMap;
	ClusterSetSubView _transformCrsClusterSets;

	Kokkos::View<TransformReactionType*> _transformReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<
		std::is_base_of_v<TransformReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = TransformReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
