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

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);
};

class FeDissociationReaction :
	public DissociationReaction<FeReactionNetwork, FeDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<FeReactionNetwork, FeDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
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
