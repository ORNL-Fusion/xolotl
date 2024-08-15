#pragma once

#include <xolotl/core/network/NEReaction.h>
#include <xolotl/core/network/NETraits.h>
#include <xolotl/core/network/ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class NEReactionGenerator;

class NEClusterUpdater;
} // namespace detail

class NEReactionNetwork : public ReactionNetwork<NEReactionNetwork>
{
	friend class ReactionNetwork<NEReactionNetwork>;

public:
	using Superclass = ReactionNetwork<NEReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	IndexType
	checkLargestClusterId();

	void
	initializeExtraClusterData(const options::IOptions& options);

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

	detail::NEReactionGenerator
	getReactionGenerator() const noexcept;

	void
	readClusters(const std::string filename = "reactionRates.txt");

	void
	readReactions(
		double temperature, const std::string filename = "reactionRates.txt");
};

namespace detail
{
class NEReactionGenerator :
	public ReactionGenerator<NEReactionNetwork, NEReactionGenerator>
{
	friend class ReactionGeneratorBase<NEReactionNetwork, NEReactionGenerator>;

public:
	using NetworkType = NEReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using Superclass =
		ReactionGenerator<NEReactionNetwork, NEReactionGenerator>;

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

class NEClusterUpdater
{
public:
	using NetworkType = NEReactionNetwork;
	using ClusterData = typename NetworkType::ClusterData;
	using IndexType = typename NetworkType::IndexType;

	KOKKOS_INLINE_FUNCTION
	void
	updateDiffusionCoefficient(const ClusterData& data, IndexType clusterId,
		IndexType gridIndex) const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/NEClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/NEReactionNetwork.tpp>
#endif
