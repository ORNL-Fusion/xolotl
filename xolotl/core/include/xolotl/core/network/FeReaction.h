#pragma once

#include <xolotl/core/network/FeTraits.h>
#include <xolotl/core/network/SinkReaction.h>

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
	computeBindingEnergy();
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
} // namespace network
} // namespace core
} // namespace xolotl
