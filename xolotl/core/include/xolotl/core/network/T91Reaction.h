#pragma once

#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/T91Traits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class T91ReactionNetwork;

class T91ProductionReaction :
	public ProductionReaction<T91ReactionNetwork, T91ProductionReaction>
{
	friend class Reaction<T91ReactionNetwork, T91ProductionReaction>;

public:
	using Superclass =
		ProductionReaction<T91ReactionNetwork, T91ProductionReaction>;

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

class T91DissociationReaction :
	public DissociationReaction<T91ReactionNetwork, T91DissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<T91ReactionNetwork, T91DissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
};

class T91SinkReaction : public SinkReaction<T91ReactionNetwork, T91SinkReaction>
{
public:
	using Superclass = SinkReaction<T91ReactionNetwork, T91SinkReaction>;

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
