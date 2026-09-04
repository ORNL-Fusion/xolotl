#pragma once

#include <xolotl/core/network/FeTraits.h>
#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/TrapReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class FeReactionNetwork;

class FeProductionReaction :
	public ProductionReaction<FeReactionNetwork, FeProductionReaction>
{
public:
	using Superclass =
		ProductionReaction<FeReactionNetwork, FeProductionReaction>;

	using Superclass::Superclass;
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
	FeProductionReaction(ReactionDataRef reactionData,
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
		// Check if the single size is involved
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
			// Reactants
			if (this->_reactants[i] < numClusters) {
				this->copyMomentIds(
					this->_reactants[i], this->_reactantMomentIds[i]);
			}
			else {
				// Bubble
					this->_reactantMomentIds[i][0] =
						this->_clusterData->heAvId(); //Placeholder of bubble
					this->_reactantMomentIds[i][1] =
						this->_clusterData->heAvId(); //for He
					this->_reactantMomentIds[i][2] =
                                                this->_clusterData->voidAvId(); //for V
			}
			// Products
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
					// Bubble
					this->_productMomentIds[i][0] =
						this->_clusterData->heAvId(); //Placeholder of bubble
					this->_productMomentIds[i][1] =
						this->_clusterData->heAvId(); // for He
					this->_productMomentIds[i][2] =
                                                this->_clusterData->voidAvId(); // for V
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
	FeProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		FeProductionReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2,
			clusterSet.cluster3)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
        void
        computeCoefficients();

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

class FeDissociationReaction :
	public DissociationReaction<FeReactionNetwork, FeDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<FeReactionNetwork, FeDissociationReaction>;

	using Superclass::Superclass;
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
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
        FeDissociationReaction(ReactionDataRef reactionData,
                const ClusterData& clusterData, IndexType reactionId,
                IndexType cluster0, IndexType cluster1, IndexType cluster2)
        {

                this->_clusterData = &clusterData;
                this->_reactionId = reactionId;
                this->_rate = reactionData.getRates(reactionId);
                this->_widths = reactionData.getWidths(reactionId);
                this->_coefs = reactionData.getCoefficients(reactionId);

                this->_reactant = cluster0;
                this->_products = {cluster1, cluster2};

                auto numClusters = clusterData.numClusters;
                // Check if the single size is involved
                if (cluster0 >= numClusters)
                        isLargeBubbleReaction = true;
                if (cluster1 != Superclass::invalidIndex and cluster1 >= numClusters)
                        isLargeBubbleReaction = true;
                if (cluster2 != Superclass::invalidIndex and cluster2 >= numClusters)
                        isLargeBubbleReaction = true;

                // static
                const auto dummyRegion = Region(Composition{});

                        // Reactant
                        if (this->_reactant < numClusters) {
                                this->copyMomentIds(
                                        this->_reactant, this->_reactantMomentIds);
                        }
                        else {
                                // Bubble
                                this->_reactantMomentIds[0] =
                                        this->_clusterData->heAvId(); // placeholeder of bubble
				this->_reactantMomentIds[1] =
                                        this->_clusterData->heAvId(); // for He
                                this->_reactantMomentIds[2] =
                                        this->_clusterData->voidAvId(); // for V
                        }
                for (auto i : {0, 1}) {
			// Products
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
					// Bubble
					this->_productMomentIds[i][0] =
						this->_clusterData->heAvId(); // placeholeder of bubble
					this->_productMomentIds[i][1] =
						this->_clusterData->heAvId(); // for He
					this->_productMomentIds[i][2] =
						this->_clusterData->voidAvId(); // for V
				}
			}
                }

		const auto& cl1Reg = (this->_reactant < numClusters) ?
                        this->_clusterData->getCluster(this->_reactant).getRegion() :
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

                this->_reactantVolume = cl1Reg.volume();
                this->_productVolumes = {pr1Reg.volume(), pr2Reg.volume()};

                this->initialize();
        }

	KOKKOS_INLINE_FUNCTION
	FeDissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		FeDissociationReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);

	KOKKOS_INLINE_FUNCTION
        void
        computeCoefficients();

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

class FeSinkReaction : public SinkReaction<FeReactionNetwork, FeSinkReaction>
{
public:
	using Superclass = SinkReaction<FeReactionNetwork, FeSinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();
};

class FeTrapReaction : public TrapReaction<FeReactionNetwork, FeTrapReaction>
{
public:
	using Superclass = TrapReaction<FeReactionNetwork, FeTrapReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	IndexType
	getId();

	KOKKOS_INLINE_FUNCTION
	double
	getEnergy();

	KOKKOS_INLINE_FUNCTION
	double
	getStrength();

	KOKKOS_INLINE_FUNCTION
	double
	getFrequency();

	KOKKOS_INLINE_FUNCTION
	double
	getDensity();
};
} // namespace network
} // namespace core
} // namespace xolotl
