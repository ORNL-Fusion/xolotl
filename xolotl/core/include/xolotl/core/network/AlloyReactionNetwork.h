#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/AlloyReaction.h>
#include <xolotl/core/network/AlloyTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class AlloyReactionGenerator;
} // namespace detail

class AlloyReactionNetwork : public ReactionNetwork<AlloyReactionNetwork>
{
	friend class ReactionNetwork<AlloyReactionNetwork>;

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

	IndexType
	getLargestClusterId()
	{
		return largestClusterId;
	}

	void
	setConstantRates(RatesView rates, IndexType gridIndex) override;

	void
	setConstantConnectivities(ConnectivitiesPair conns) override;

	void
	setConstantRateEntries() override;

	std::string
	getMonitorOutputFileName() const override
	{
		return "Alloy.dat";
	}

	void
	initializeExtraClusterData(const options::IOptions& options);

	void
	initializeExtraDOFs(const options::IOptions& options);

	void
	computeFluxesPreProcess(ConcentrationsView concentrations,
		FluxesView fluxes, IndexType gridIndex, double surfaceDepth,
		double spacing);

	void
	computePartialsPreProcess(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
		double spacing);

	double
	computeBubbleRadius(
		double amount, double latticeParameter, IndexType typeSwitch = 0)
	{
		// Find the edge of the phase space
		const auto& largestReg = this->getCluster(largestClusterId).getRegion();
		Composition hiLargest = largestReg.getUpperLimitPoint();
		double largestSize = hiLargest[Species::He] + hiLargest[Species::V] +
			hiLargest[Species::PerfectV] + hiLargest[Species::FaultedV] +
			hiLargest[Species::PerfectI] + hiLargest[Species::FaultedI] -
			5; // Don't know which one was saved

		// Get the minimum amount for a valid value
		amount = util::max(amount, largestSize);

		// Prefactor
		const double prefactor =
			0.25 * latticeParameter * latticeParameter / ::xolotl::core::pi;

		// Compute the radius
		switch (typeSwitch) {
		// Bubble
		case 0:
			return cbrt(0.75 * prefactor * latticeParameter * amount);
		// Perfect V
		case 1:
			return sqrt((amount * prefactor) / ::xolotl::core::perfectBurgers);
		// Faulted V
		case 2:
			return sqrt((amount * prefactor) / ::xolotl::core::faultedBurgers);
		// Perfect I
		case 3:
			return sqrt((amount * prefactor) / ::xolotl::core::perfectBurgers);
		// Faulted I
		case 4:
			return sqrt((amount * prefactor) / ::xolotl::core::faultedBurgers);
		}
		return 0.0;
	}

	std::string
	getMonitorDataHeaderString() const override;

	void
	addMonitorDataValues(Kokkos::View<const double*> conc, double fac,
		std::vector<double>& totalVals) override;

	std::size_t
	getMonitorDataLineSize() const override
	{
		return getSpeciesListSize() * 4;
	}

	void
	writeMonitorDataLine(
		const std::vector<double>& localData, double time) override;

public:
	IndexType largestClusterId;

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

	AlloyReactionGenerator(const AlloyReactionNetwork& network) :
		Superclass(network),
		largestClusterId(network.largestClusterId)
	{
	}

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

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSingleSizeReactions(IndexType i, IndexType j, TTag tag) const;

private:
	ReactionCollection<Network>
	getReactionCollection() const;

	IndexType largestClusterId;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/AlloyClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/AlloyReactionNetwork.tpp>
#endif
