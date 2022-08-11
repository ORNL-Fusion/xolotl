#pragma once

#include <xolotl/core/network/FeReaction.h>
#include <xolotl/core/network/FeTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class FeReactionGenerator;
}

class FeReactionNetwork : public ReactionNetwork<FeReactionNetwork>
{
	friend class ReactionNetwork<FeReactionNetwork>;
	friend class detail::ReactionNetworkWorker<FeReactionNetwork>;

public:
	using Superclass = ReactionNetwork<FeReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using AmountType = typename Superclass::AmountType;
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
		// 2 atoms per cell
		return 0.5 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::FeReactionGenerator
	getReactionGenerator() const noexcept;

	void
	readClusters(const std::string filename)
	{
		return;
	}

	void
	readReactions(double temperature, const std::string filename)
	{
		return;
	}
};

namespace detail
{
class FeReactionGenerator :
	public ReactionGenerator<FeReactionNetwork, FeReactionGenerator>
{
	friend class ReactionGeneratorBase<FeReactionNetwork, FeReactionGenerator>;

public:
	using NetworkType = FeReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;

	using Superclass =
		ReactionGenerator<FeReactionNetwork, FeReactionGenerator>;

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
	ReactionCollection<NetworkType>
	getReactionCollection() const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/FeClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/FeReactionNetwork.tpp>
#endif
