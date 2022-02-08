#pragma once

#include <xolotl/core/network/AlloyReaction.h>
#include <xolotl/core/network/AlloyTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class AlloyReactionGenerator;

class AlloyClusterUpdater;
} // namespace detail

class AlloyReactionNetwork : public ReactionNetwork<AlloyReactionNetwork>
{
	friend class ReactionNetwork<AlloyReactionNetwork>;
	friend class detail::ReactionNetworkWorker<AlloyReactionNetwork>;

public:
	using Superclass = ReactionNetwork<AlloyReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	IndexType
	checkLargestClusterId();

private:
	double
	checkLatticeParameter(double latticeParameter);

	double
	computeAtomicVolume(double latticeParameter)
	{
		// 4 atoms per cell
		return 0.25 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::AlloyReactionGenerator
	getReactionGenerator() const noexcept;

	void
	readReactions(double temperature, const std::string filename)
	{
		return;
	}
};

namespace detail
{
class AlloyReactionGenerator :
	public ReactionGenerator<AlloyReactionNetwork, AlloyReactionGenerator>
{
	friend class ReactionGeneratorBase<AlloyReactionNetwork,
		AlloyReactionGenerator>;

public:
	using Network = AlloyReactionNetwork;
	using Subpaving = typename Network::Subpaving;
	using Superclass =
		ReactionGenerator<AlloyReactionNetwork, AlloyReactionGenerator>;

	using Superclass::Superclass;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSinks(IndexType i, TTag tag) const;

private:
	ReactionCollection<Network>
	getReactionCollection() const;
};

class AlloyClusterUpdater
{
public:
	using Network = AlloyReactionNetwork;
	using ClusterData = typename Network::ClusterData;
	using IndexType = typename Network::IndexType;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/AlloyClusterGenerator.h>
