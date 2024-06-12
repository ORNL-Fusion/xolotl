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
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/AlloyClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/AlloyReactionNetwork.tpp>
#endif
