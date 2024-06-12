#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/core/network/ZrReaction.h>
#include <xolotl/core/network/ZrTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class ZrReactionGenerator;

class ZrClusterUpdater;
} // namespace detail

class ZrReactionNetwork : public ReactionNetwork<ZrReactionNetwork>
{
	friend class ReactionNetwork<ZrReactionNetwork>;
	friend class detail::ReactionNetworkWorker<ZrReactionNetwork>;

public:
	using Superclass = ReactionNetwork<ZrReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using RateVector = typename Superclass::RateVector;
	using ConnectivitiesPair = typename Superclass::ConnectivitiesPair;
	using RatesView = typename Superclass::RatesView;

	using Superclass::Superclass;

	IndexType
	checkLargestClusterId();

	void
	setConstantRates(RatesView rates, IndexType gridIndex) override;

	void
	setConstantConnectivities(ConnectivitiesPair conns) override;

	void
	setConstantRateEntries() override;

	void
	initializeExtraClusterData(const options::IOptions& options);

	void
	setGridSize(IndexType gridSize) override;

	std::string
	getMonitorOutputFileName() const override
	{
		return "AlphaZr.dat";
	}

	std::string
	getMonitorDataHeaderString() const override;

	void
	addMonitorDataValues(Kokkos::View<const double*> conc, double fac,
		std::vector<double>& totalVals) override;

	std::size_t
	getMonitorDataLineSize() const override
	{
		return getSpeciesListSize() * 6;
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
		// sqrt(3) / 4 * a^2 * c
		// with c the other lattice parameter
		return 0.0234; // nm^3
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::ZrReactionGenerator
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

	void
	defineReactions(Connectivity& connectivity);
};

namespace detail
{
class ZrReactionGenerator :
	public ReactionGenerator<ZrReactionNetwork, ZrReactionGenerator>
{
	friend class ReactionGeneratorBase<ZrReactionNetwork, ZrReactionGenerator>;

public:
	using Network = ZrReactionNetwork;
	using Subpaving = typename Network::Subpaving;
	using Superclass =
		ReactionGenerator<ZrReactionNetwork, ZrReactionGenerator>;

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

class ZrClusterUpdater
{
public:
	using Network = ZrReactionNetwork;
	using ClusterData = typename Network::ClusterData;
	using IndexType = typename Network::IndexType;

	KOKKOS_INLINE_FUNCTION
	void
	updateDiffusionCoefficient(const ClusterData& data, IndexType clusterId,
		IndexType gridIndex) const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/ZrClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/ZrReactionNetwork.tpp>
#endif
