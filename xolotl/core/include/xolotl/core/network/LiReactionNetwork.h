#pragma once

#include <xolotl/core/network/LiReaction.h>
#include <xolotl/core/network/LiTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class LiReactionGenerator;
}

class LiReactionNetwork : public ReactionNetwork<LiReactionNetwork>
{
	friend class ReactionNetwork<LiReactionNetwork>;

public:
	using Superclass = ReactionNetwork<LiReactionNetwork>;
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
		// TODO
		return 0.5 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::LiReactionGenerator
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
class LiReactionGenerator :
	public ReactionGenerator<LiReactionNetwork, LiReactionGenerator>
{
	friend class ReactionGeneratorBase<LiReactionNetwork, LiReactionGenerator>;

public:
	using NetworkType = LiReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;

	using Superclass =
		ReactionGenerator<LiReactionNetwork, LiReactionGenerator>;

	using Superclass::Superclass;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/LiClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/LiReactionNetwork.tpp>
#endif
