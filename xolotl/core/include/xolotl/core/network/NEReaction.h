#pragma once

#include <xolotl/core/network/NETraits.h>
#include <xolotl/core/network/NucleationReaction.h>
#include <xolotl/core/network/ReSolutionReaction.h>
#include <xolotl/core/network/SinkReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class NEReactionNetwork;

class NEProductionReaction :
	public ProductionReaction<NEReactionNetwork, NEProductionReaction>
{
	friend class Reaction<NEReactionNetwork, NEProductionReaction>;

public:
	using Superclass =
		ProductionReaction<NEReactionNetwork, NEProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	NEProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1,
		IndexType cluster2 = invalidIndex, IndexType cluster3 = invalidIndex) :
		ProductionReaction(reactionData, clusterData, reactionId, cluster0,
			cluster1, cluster2, cluster3)
	{
		if (this->_clusterData->extraData.fileClusterMap.exists(
				_reactants[0]) and
			this->_clusterData->extraData.fileClusterMap.exists(
				_reactants[1])) {
			auto clusterMap = this->_clusterData->extraData.fileClusterMap;
			this->_deltaG0 = reactionData.reactionEnergies(
				clusterMap.value_at(clusterMap.find(_reactants[0])),
				clusterMap.value_at(clusterMap.find(_reactants[1])));
		}
	}

	KOKKOS_INLINE_FUNCTION
	NEProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		NEProductionReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2,
			clusterSet.cluster3)
	{
	}

private:
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

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);
};

class NEDissociationReaction :
	public DissociationReaction<NEReactionNetwork, NEDissociationReaction>
{
	friend class Reaction<NEReactionNetwork, NEDissociationReaction>;

public:
	using Superclass =
		DissociationReaction<NEReactionNetwork, NEDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	NEDissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1, IndexType cluster2) :
		DissociationReaction(
			reactionData, clusterData, reactionId, cluster0, cluster1, cluster2)
	{
		if (this->_clusterData->extraData.fileClusterMap.exists(
				_products[0]) and
			this->_clusterData->extraData.fileClusterMap.exists(_products[1])) {
			auto clusterMap = this->_clusterData->extraData.fileClusterMap;
			this->_deltaG0 = reactionData.reactionEnergies(
				clusterMap.value_at(clusterMap.find(_products[0])),
				clusterMap.value_at(clusterMap.find(_products[1])));
		}
	}

	KOKKOS_INLINE_FUNCTION
	NEDissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		NEDissociationReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy();

private:
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

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);
};

class NEReSolutionReaction :
	public ReSolutionReaction<NEReactionNetwork, NEReSolutionReaction>
{
public:
	using Superclass =
		ReSolutionReaction<NEReactionNetwork, NEReSolutionReaction>;

	using Superclass::Superclass;
};

class NENucleationReaction :
	public NucleationReaction<NEReactionNetwork, NENucleationReaction>
{
public:
	using Superclass =
		NucleationReaction<NEReactionNetwork, NENucleationReaction>;

	using Superclass::Superclass;
};

class NESinkReaction : public SinkReaction<NEReactionNetwork, NESinkReaction>
{
	friend class Reaction<NEReactionNetwork, NESinkReaction>;

public:
	using Superclass = SinkReaction<NEReactionNetwork, NESinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	NESinkReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, IndexType cluster0) :
		SinkReaction(reactionData, clusterData, reactionId, cluster0)
	{
		if (this->_clusterData->extraData.fileClusterMap.exists(_reactant)) {
			auto clusterMap = this->_clusterData->extraData.fileClusterMap;
			this->_deltaG0 = reactionData.reactionEnergies(
				clusterMap.value_at(clusterMap.find(_reactant)),
				clusterMap.capacity());
		}
	}

	KOKKOS_INLINE_FUNCTION
	NESinkReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet) :
		NESinkReaction(
			reactionData, clusterData, reactionId, clusterSet.cluster0)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();
};
} // namespace network
} // namespace core
} // namespace xolotl
