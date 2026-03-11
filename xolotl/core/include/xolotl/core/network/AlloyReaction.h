#pragma once

#include <xolotl/core/network/AlloyTraits.h>
#include <xolotl/core/network/ConstantReaction.h>
#include <xolotl/core/network/NucleationReaction.h>
#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/TransformReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class AlloyReactionNetwork;

class AlloyProductionReaction :
	public ProductionReaction<AlloyReactionNetwork, AlloyProductionReaction>
{
public:
	using Superclass =
		ProductionReaction<AlloyReactionNetwork, AlloyProductionReaction>;

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
	AlloyProductionReaction(ReactionDataRef reactionData,
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
				auto shift = (this->_reactants[i] - numClusters) / 2;
				switch (shift) {
				// Bubble
				case 0:
					this->_reactantMomentIds[i][0] =
						this->_clusterData->voidAvId();
					this->_reactantMomentIds[i][1] =
						this->_clusterData->heAvId();
					for (IndexType j = 2; j < Superclass::nMomentIds; ++j) {
						this->_reactantMomentIds[i][j] =
							Superclass::invalidIndex;
					}
					break;
				// Perfect V
				case 1:
					this->_reactantMomentIds[i][0] =
						this->_clusterData->perfVAvId();
					for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
						this->_reactantMomentIds[i][j] =
							Superclass::invalidIndex;
					}
					break;
				// Faulted V
				case 2:
					this->_reactantMomentIds[i][0] =
						this->_clusterData->faulVAvId();
					for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
						this->_reactantMomentIds[i][j] =
							Superclass::invalidIndex;
					}
					break;
				// Perfect I
				case 3:
					this->_reactantMomentIds[i][0] =
						this->_clusterData->perfIAvId();
					for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
						this->_reactantMomentIds[i][j] =
							Superclass::invalidIndex;
					}
					break;
				// Faulted I
				case 4:
					this->_reactantMomentIds[i][0] =
						this->_clusterData->faulIAvId();
					for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
						this->_reactantMomentIds[i][j] =
							Superclass::invalidIndex;
					}
					break;
				}
			}
			// Products
			if (this->_products[i] == Superclass::invalidIndex) {
				for (IndexType j = 0; j < Superclass::nMomentIds; ++j) {
					this->_productMomentIds[i][j] = Superclass::invalidIndex;
				}
			}
			else {
				if (this->_products[i] < numClusters) {
					this->copyMomentIds(
						this->_products[i], this->_productMomentIds[i]);
				}
				else {
					auto shift = (this->_products[i] - numClusters) / 2;
					switch (shift) {
					// Bubble
					case 0:
						this->_productMomentIds[i][0] =
							this->_clusterData->voidAvId();
						this->_productMomentIds[i][1] =
							this->_clusterData->heAvId();
						for (IndexType j = 2; j < Superclass::nMomentIds; ++j) {
							this->_productMomentIds[i][j] =
								Superclass::invalidIndex;
						}
						break;
					// Perfect V
					case 1:
						this->_productMomentIds[i][0] =
							this->_clusterData->perfVAvId();
						for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
							this->_productMomentIds[i][j] =
								Superclass::invalidIndex;
						}
						break;
					// Faulted V
					case 2:
						this->_productMomentIds[i][0] =
							this->_clusterData->faulVAvId();
						for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
							this->_productMomentIds[i][j] =
								Superclass::invalidIndex;
						}
						break;
					// Perfect I
					case 3:
						this->_productMomentIds[i][0] =
							this->_clusterData->perfIAvId();
						for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
							this->_productMomentIds[i][j] =
								Superclass::invalidIndex;
						}
						break;
					// Faulted I
					case 4:
						this->_productMomentIds[i][0] =
							this->_clusterData->faulIAvId();
						for (IndexType j = 1; j < Superclass::nMomentIds; ++j) {
							this->_productMomentIds[i][j] =
								Superclass::invalidIndex;
						}
						break;
					}
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
	AlloyProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		AlloyProductionReaction(reactionData, clusterData, reactionId,
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

class AlloyDissociationReaction :
	public DissociationReaction<AlloyReactionNetwork, AlloyDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<AlloyReactionNetwork, AlloyDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
};

class AlloySinkReaction :
	public SinkReaction<AlloyReactionNetwork, AlloySinkReaction>
{
public:
	using Superclass = SinkReaction<AlloyReactionNetwork, AlloySinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();
};

class AlloyTransformReaction :
	public TransformReaction<AlloyReactionNetwork, AlloyTransformReaction>
{
public:
	using Superclass =
		TransformReaction<AlloyReactionNetwork, AlloyTransformReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	AlloyTransformReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1)
	{
		this->_clusterData = &clusterData;
		this->_reactionId = reactionId;
		this->_rate = reactionData.getRates(reactionId);
		this->_widths = reactionData.getWidths(reactionId);
		this->_coefs = reactionData.getCoefficients(reactionId);

		this->_reactant = cluster0;
		this->_product = cluster1;

		auto numClusters = clusterData.numClusters;
		// Check if the single size is involved
		if (cluster0 >= numClusters)
			isLargeBubbleReaction = true;

		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	AlloyTransformReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		AlloyTransformReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getSize();

	KOKKOS_INLINE_FUNCTION
	double
	getExponent();

	KOKKOS_INLINE_FUNCTION
	double
	getBarrier();

	KOKKOS_INLINE_FUNCTION
	void
	computeConnectivity(const Connectivity& connectivity)
	{
		Superclass::computeConnectivity(connectivity);

		if (isLargeBubbleReaction) {
			this->addConnectivity(
				this->_product + 1, this->_reactant, connectivity);
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	mapJacobianEntries(Connectivity connectivity)
	{
		Superclass::mapJacobianEntries(connectivity);

		if (isLargeBubbleReaction) {
			_connEntries[1][1][0][0] = connectivity(_product + 1, _reactant);
		}
	}

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

class AlloyNucleationReaction :
	public NucleationReaction<AlloyReactionNetwork, AlloyNucleationReaction>
{
public:
	using Superclass =
		NucleationReaction<AlloyReactionNetwork, AlloyNucleationReaction>;

	using Superclass::Superclass;
};

class AlloyConstantReaction :
	public ConstantReaction<AlloyReactionNetwork, AlloyConstantReaction>
{
public:
	using Superclass =
		ConstantReaction<AlloyReactionNetwork, AlloyConstantReaction>;

	using Superclass::Superclass;
};
} // namespace network
} // namespace core
} // namespace xolotl
