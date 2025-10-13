#pragma once

#include <xolotl/core/network/RPVTraits.h>
#include <xolotl/core/network/SinkReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class RPVReactionNetwork;

class RPVProductionReaction :
	public ProductionReaction<RPVReactionNetwork, RPVProductionReaction>
{
	friend class Reaction<RPVReactionNetwork, RPVProductionReaction>;

public:
	using Superclass =
		ProductionReaction<RPVReactionNetwork, RPVProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

private:
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

class RPVDissociationReaction :
	public DissociationReaction<RPVReactionNetwork, RPVDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<RPVReactionNetwork, RPVDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
};

class RPVSinkReaction : public SinkReaction<RPVReactionNetwork, RPVSinkReaction>
{
public:
	using Superclass = SinkReaction<RPVReactionNetwork, RPVSinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);
};
} // namespace network
} // namespace core
} // namespace xolotl
