#pragma once

#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/PSIReaction.h>
#include <xolotl/core/network/PSITraits.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/core/network/detail/ReactionGenerator.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TSpeciesEnum>
class PSIReactionGenerator;
}

template <typename TSpeciesEnum>
struct ReactionNetworkInterface<PSIReactionNetwork<TSpeciesEnum>>
{
	using Type = IPSIReactionNetwork;
};

template <typename TSpeciesEnum>
class PSIReactionNetwork :
	public ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>
{
	friend class ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
	friend class detail::ReactionNetworkWorker<
		PSIReactionNetwork<TSpeciesEnum>>;

public:
	using Superclass = ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using AmountType = typename Superclass::AmountType;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	SpeciesId
	getHeliumSpeciesId() const override
	{
		return SpeciesId{Species::He, this->getNumberOfSpecies()};
	}

	SpeciesId
	getInterstitialSpeciesId() const override
	{
		return SpeciesId{Species::I, this->getNumberOfSpecies()};
	}

	SpeciesId
	getVacancySpeciesId() const override
	{
		return SpeciesId{Species::V, this->getNumberOfSpecies()};
	}

	double
	getTotalTrappedHeliumConcentration(
		ConcentrationsView concs, AmountType minSize = 0) override
	{
		return this->getTotalTrappedAtomConcentration(
			concs, Species::He, minSize);
	}

	void
	updateBurstingConcs(double* gridPointSolution, double factor,
		std::vector<double>& nBurst) override;

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

	detail::PSIReactionGenerator<Species>
	getReactionGenerator() const noexcept
	{
		return detail::PSIReactionGenerator<Species>{*this};
	}
};

namespace detail
{
template <typename TSpeciesEnum>
class PSIReactionGenerator :
	public ReactionGenerator<PSIReactionNetwork<TSpeciesEnum>,
		PSIReactionGenerator<TSpeciesEnum>>
{
	friend class ReactionGeneratorBase<PSIReactionNetwork<TSpeciesEnum>,
		PSIReactionGenerator<TSpeciesEnum>>;

public:
	using NetworkType = PSIReactionNetwork<TSpeciesEnum>;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;

	using Superclass = ReactionGenerator<PSIReactionNetwork<TSpeciesEnum>,
		PSIReactionGenerator<TSpeciesEnum>>;

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

#include <xolotl/core/network/PSIClusterGenerator.h>
