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
template <typename TBase>
class NucleationReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using NucleationReactionType =
		typename NetworkTraits::NucleationReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	NucleationReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addNucleationReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addNucleationReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<NucleationReactionType*>
	getNucleationReactions() const
	{
		return _nucleationReactions;
	}

	IndexType
	getNumberOfNucleationReactions() const
	{
		return _nucleationReactions.size();
	}

private:
	IndexView _clusterNucleationReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numNucleationReactions{};

	IndexView _nucleationCrsRowMap;
	ClusterSetSubView _nucleationCrsClusterSets;

	Kokkos::View<NucleationReactionType*> _nucleationReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<std::is_base_of<NucleationReaction<TNetwork, TReaction>,
		TReaction>::value>>
{
	using Type = NucleationReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
