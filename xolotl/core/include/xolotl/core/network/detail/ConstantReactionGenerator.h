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
 * constant reactions depending on the subpaving.
 *
 * @tparam TBase The templated base class.
 */
template <typename TBase>
class ConstantReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using ConstantReactionType = typename NetworkTraits::ConstantReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	ConstantReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addConstantReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addConstantReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<ConstantReactionType*>
	getConstantReactions() const
	{
		return _constantReactions;
	}

	IndexType
	getNumberOfConstantReactions() const
	{
		return _constantReactions.size();
	}

private:
	IndexView _clusterConstantReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numConstantReactions{};

	IndexView _constantCrsRowMap;
	ClusterSetSubView _constantCrsClusterSets;

	Kokkos::View<ConstantReactionType*> _constantReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<
		std::is_base_of_v<ConstantReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = ConstantReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
