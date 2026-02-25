#pragma once

#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/VTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class VReactionNetwork;

class VProductionReaction :
	public ProductionReaction<VReactionNetwork, VProductionReaction>
{
public:
	using Superclass =
		ProductionReaction<VReactionNetwork, VProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);
};

class VDissociationReaction :
	public DissociationReaction<VReactionNetwork, VDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<VReactionNetwork, VDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
};

class VSinkReaction : public SinkReaction<VReactionNetwork, VSinkReaction>
{
public:
	using Superclass = SinkReaction<VReactionNetwork, VSinkReaction>;

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
