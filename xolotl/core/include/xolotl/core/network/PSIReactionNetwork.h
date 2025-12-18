#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/PSIReaction.h>
#include <xolotl/core/network/PSITraits.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/core/network/detail/ReactionGenerator.h>
#include <xolotl/core/network/detail/TrapMutationHandler.h>
#include <xolotl/util/MathUtils.h>

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

public:
	using Superclass = ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
	using Subpaving = typename Superclass::Subpaving;
	using SubdivisionRatio = typename Superclass::SubdivisionRatio;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using AmountType = typename Superclass::AmountType;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	PSIReactionNetwork(const Subpaving& subpaving, IndexType gridSize,
		const options::IOptions& options);

	PSIReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		const std::vector<SubdivisionRatio>& subdivisionRatios,
		IndexType gridSize, const options::IOptions& options);

	PSIReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
		IndexType gridSize, const options::IOptions& options);

	SpeciesId
	getHeliumSpeciesId() const override
	{
		return SpeciesId{Species::He, Superclass::getNumberOfSpecies()};
	}

	SpeciesId
	getVacancySpeciesId() const override
	{
		return SpeciesId{Species::V, Superclass::getNumberOfSpecies()};
	}

	SpeciesId
	getInterstitialSpeciesId() const override
	{
		return SpeciesId{Species::I, Superclass::getNumberOfSpecies()};
	}

	bool
	hasDeuterium() const noexcept override
	{
		return psi::hasDeuterium<Species>;
	}

	bool
	hasTritium() const noexcept override
	{
		return psi::hasTritium<Species>;
	}

	void
	initializeExtraClusterData(const options::IOptions& options);

	void
	updateExtraClusterData(const std::vector<double>& gridTemps,
		const std::vector<double>& gridDepths);

	void
	initializeExtraDOFs(const options::IOptions& options);

	void
	selectTrapMutationReactions(double surfaceDepth, double spacing);

	void
	computeFluxesPreProcess(ConcentrationsView concentrations,
		FluxesView fluxes, IndexType gridIndex, double surfaceDepth,
		double spacing);

	void
	computePartialsPreProcess(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
		double spacing);

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

	double
	computeBubbleRadius(double amount, double latticeParameter)
	{
		// Find the edge of the phase space
		const auto& largestReg = this->getCluster(largestClusterId).getRegion();
		Composition hiLargest = largestReg.getUpperLimitPoint();
		double largestSize = hiLargest[Species::V] - 1;

		// Get the minimum amount for a valid value
		amount = util::max(amount, largestSize);

		return (sqrt(3.0) / 4.0) * latticeParameter +
			pow((3.0 * pow(latticeParameter, 3.0) * amount) /
					(8.0 * ::xolotl::core::pi),
				(1.0 / 3.0)) -
			pow((3.0 * pow(latticeParameter, 3.0)) / (8.0 * ::xolotl::core::pi),
				(1.0 / 3.0));
	}

	IndexType
	checkLargestClusterId();

	IndexType
	getLargestClusterId()
	{
		return largestClusterId;
	}

	void
	updateReactionRates(double time = 0.0);

	void
	updateTrapMutationDisappearingRate(double totalTrappedHeliumConc) override;

	void
	updateDesorptionLeftSideRate(
		ConcentrationsView concentrations, IndexType gridIndex);

	// Save the ID of the largest cluster in the network
	IndexType largestClusterId;

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
	getReactionGenerator() noexcept
	{
		return detail::PSIReactionGenerator<Species>{*this};
	}

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

private:
	std::unique_ptr<detail::TrapMutationHandler> _tmHandler;
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
	using AmountType = typename NetworkType::AmountType;

	using Superclass = ReactionGenerator<PSIReactionNetwork<TSpeciesEnum>,
		PSIReactionGenerator<TSpeciesEnum>>;

	PSIReactionGenerator(PSIReactionNetwork<TSpeciesEnum>& network);

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
	addSingleSizeReactions(IndexType i, IndexType j, TTag tag) const;

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;

private:
	Kokkos::Array<Kokkos::View<AmountType*>, 7> _tmVSizes;

	bool hasHelium = false;

	// Save the ID of the largest cluster in the network
	IndexType largestClusterId;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/PSIClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/PSIReactionNetwork.tpp>
#endif
