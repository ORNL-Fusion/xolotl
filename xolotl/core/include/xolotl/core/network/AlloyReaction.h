#pragma once

#include <xolotl/core/network/AlloyTraits.h>
#include <xolotl/core/network/SinkReaction.h>

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
	computeBindingEnergy();
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
};
} // namespace network
} // namespace core
} // namespace xolotl
