#pragma once

#include <xolotl/core/network/BurstingReaction.h>
#include <xolotl/core/network/PSITraits.h>
#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/TrapMutationReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIProductionReaction :
	public ProductionReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIProductionReaction<TSpeciesEnum>>
{
public:
	using Superclass = ProductionReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIProductionReaction<TSpeciesEnum>>;

	using Superclass::Superclass;
	using Connectivity = typename Superclass::Connectivity;
	using NetworkType = typename Superclass::NetworkType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;
	using IndexType = typename Superclass::IndexType;
	using Composition = typename Superclass::Composition;
	using Region = typename Superclass::Region;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using Species = typename Superclass::Species;

	KOKKOS_INLINE_FUNCTION
	PSIProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1,
		IndexType cluster2 = Superclass::invalidIndex,
		IndexType cluster3 = Superclass::invalidIndex)
	{
		this->_clusterData = &clusterData;
		this->_reactionId = reactionId;
		this->_rate = reactionData.getRates(reactionId);
		this->_widths = reactionData.getWidths(reactionId);
		this->_coefs = reactionData.getCoefficients(reactionId);

		this->_reactants = {cluster0, cluster1};
		this->_products = {cluster2, cluster3};

		auto numClusters = clusterData.numClusters;
		// Check if the large bubble is involved
		if (cluster0 >= numClusters)
			isLargeBubbleReaction = true;
		if (cluster1 >= numClusters)
			isLargeBubbleReaction = true;
		if (cluster2 != Superclass::invalidIndex and cluster2 >= numClusters)
			isLargeBubbleReaction = true;
		if (cluster3 != Superclass::invalidIndex and cluster3 >= numClusters)
			isLargeBubbleReaction = true;

		// static
		const auto dummyRegion = Region(Composition{});

		for (auto i : {0, 1}) {
			if (this->_reactants[i] < numClusters) {
				this->copyMomentIds(
					this->_reactants[i], this->_reactantMomentIds[i]);
			}
			else {
				this->_reactantMomentIds[i][0] =
					this->_clusterData->bubbleAvHeId();
				this->_reactantMomentIds[i][1] =
					this->_clusterData->bubbleAvVId();
			}
			if (this->_products[i] < numClusters) {
				this->copyMomentIds(
					this->_products[i], this->_productMomentIds[i]);
			}
			else {
				if (this->_products[i] == Superclass::invalidIndex) {
					for (IndexType j = 0; j < Superclass::nMomentIds; ++j) {
						this->_productMomentIds[i][j] =
							Superclass::invalidIndex;
					}
				}
				else {
					this->_productMomentIds[i][0] =
						this->_clusterData->bubbleAvHeId();
					this->_productMomentIds[i][1] =
						this->_clusterData->bubbleAvVId();
				}
			}
		}

		const auto& cl1Reg = (this->_reactants[0] < numClusters) ?
			this->_clusterData->getCluster(this->_reactants[0]).getRegion() :
			dummyRegion;
		const auto& cl2Reg = (this->_reactants[1] < numClusters) ?
			this->_clusterData->getCluster(this->_reactants[1]).getRegion() :
			dummyRegion;
		const auto& pr1Reg = (this->_products[0] == Superclass::invalidIndex) ?
			dummyRegion :
			(this->_products[0] < numClusters) ?
			this->_clusterData->getCluster(this->_products[0]).getRegion() :
			dummyRegion;
		const auto& pr2Reg = (this->_products[1] == Superclass::invalidIndex) ?
			dummyRegion :
			(this->_products[1] < numClusters) ?
			this->_clusterData->getCluster(this->_products[1]).getRegion() :
			dummyRegion;

		this->_reactantVolumes = {cl1Reg.volume(), cl2Reg.volume()};
		this->_productVolumes = {pr1Reg.volume(), pr2Reg.volume()};

		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	PSIProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		PSIProductionReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2,
			clusterSet.cluster3)
	{
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients();

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

private:
	bool isLargeBubbleReaction = false;
};

template <typename TSpeciesEnum>
class PSIDissociationReaction :
	public DissociationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIDissociationReaction<TSpeciesEnum>>
{
public:
	using Superclass = DissociationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIDissociationReaction<TSpeciesEnum>>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy();
};

template <typename TSpeciesEnum>
class PSISinkReaction :
	public SinkReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSISinkReaction<TSpeciesEnum>>
{
public:
	using Superclass = SinkReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSISinkReaction<TSpeciesEnum>>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();
};

template <typename TSpeciesEnum>
class PSITrapMutationReaction :
	public TrapMutationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSITrapMutationReaction<TSpeciesEnum>>
{
public:
	using Superclass = TrapMutationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSITrapMutationReaction<TSpeciesEnum>>;
	using Superclass::Superclass;
};
template <typename TSpeciesEnum>
class PSIBurstingReaction :
	public BurstingReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIBurstingReaction<TSpeciesEnum>>
{
public:
	using Superclass = BurstingReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIBurstingReaction<TSpeciesEnum>>;
	using IndexType = typename Superclass::IndexType;
	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getAppliedRate(IndexType gridIndex) const;
};
} // namespace network
} // namespace core
} // namespace xolotl
