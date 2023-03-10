#pragma once

#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/core/network/FeCrReaction.h>
#include <xolotl/core/network/FeCrTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class FeCrReactionGenerator;

class FeCrClusterUpdater;
} // namespace detail

class FeCrReactionNetwork : public ReactionNetwork<FeCrReactionNetwork>
{
	friend class ReactionNetwork<FeCrReactionNetwork>;
	friend class detail::ReactionNetworkWorker<FeCrReactionNetwork>;

public:
	using Superclass = ReactionNetwork<FeCrReactionNetwork>;
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

	void
	computeFluxesPreProcess(ConcentrationsView concentrations,
		FluxesView fluxes, IndexType gridIndex, double surfaceDepth,
		double spacing);

	void
	computePartialsPreProcess(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
		double spacing);

	void
	updateOutgoingSinkFluxes(double* gridPointSolution,
		std::vector<double>& fluxes, IndexType gridIndex) override;

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

	detail::FeCrReactionGenerator
	getReactionGenerator() const noexcept;
};

namespace detail
{
class FeCrReactionGenerator :
	public ReactionGenerator<FeCrReactionNetwork, FeCrReactionGenerator>
{
	friend class ReactionGeneratorBase<FeCrReactionNetwork,
		FeCrReactionGenerator>;

public:
	using Network = FeCrReactionNetwork;
	using Subpaving = typename Network::Subpaving;
	using Superclass =
		ReactionGenerator<FeCrReactionNetwork, FeCrReactionGenerator>;

	using Superclass::Superclass;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSinks(IndexType i, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addTransforms(IndexType i, IndexType j, TTag tag) const;

private:
	ReactionCollection<Network>
	getReactionCollection() const;
};

class FeCrClusterUpdater
{
public:
	using Network = FeCrReactionNetwork;
	using ClusterData = typename Network::ClusterData;
	using IndexType = typename Network::IndexType;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/FeCrClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/FeCrReactionNetwork.tpp>
#endif
