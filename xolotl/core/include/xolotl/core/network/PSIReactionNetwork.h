#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/core/network/PSIReaction.h>
#include <xolotl/core/network/PSITraits.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/core/network/detail/ReactionGenerator.h>
#include <xolotl/core/network/detail/TrapMutationHandler.h>

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
	using SubdivisionRatio = typename Superclass::SubdivisionRatio;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using AmountType = typename Superclass::AmountType;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using Connectivity = typename Superclass::Connectivity;

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
	initializeExtraDOFs(const options::IOptions& options);

	void
	updateExtraClusterData(const std::vector<double>& gridTemps,
		const std::vector<double>& gridDepths);

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
		auto toReturn =
			this->getTotalTrappedAtomConcentration(concs, Species::He, minSize);

		if (this->_enableLargeBubble) {
			// Get the large bubble contribution
			auto dof = concs.size();
			auto dHe = Kokkos::subview(concs, std::make_pair(dof - 2, dof - 1));
			auto hHe = create_mirror_view(dHe);
			deep_copy(hHe, dHe);

			toReturn += hHe(0);
		}

		return toReturn;
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

private:
	std::unique_ptr<detail::TrapMutationHandler> _tmHandler;

protected:
	double
	computeBubbleRadius(double vAmount, double latticeParameter)
	{
		double aCube = latticeParameter * latticeParameter * latticeParameter;

		return (sqrt(3.0) / 4.0) * latticeParameter +
			pow((3.0 * aCube * vAmount) / (8.0 * ::xolotl::core::pi),
				(1.0 / 3.0)) -
			pow((3.0 * aCube) / (8.0 * ::xolotl::core::pi), (1.0 / 3.0));
	}

public:
	IndexType largestClusterId;
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

	PSIReactionGenerator(const PSIReactionNetwork<TSpeciesEnum>& network);

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
	addLargeBubbleReactions(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addBurstings(IndexType i, TTag tag) const;

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;

private:
	Kokkos::Array<Kokkos::View<AmountType*>, 7> _tmVSizes;

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
