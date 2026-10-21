#pragma once

#include <xolotl/core/network/LiTraits.h>
#include <xolotl/core/network/SinkReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class LiReactionNetwork;

class LiProductionReaction :
	public ProductionReaction<LiReactionNetwork, LiProductionReaction>
{
public:
	using Superclass =
		ProductionReaction<LiReactionNetwork, LiProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);
};

class LiDissociationReaction :
	public DissociationReaction<LiReactionNetwork, LiDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<LiReactionNetwork, LiDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
};
} // namespace network
} // namespace core
} // namespace xolotl
