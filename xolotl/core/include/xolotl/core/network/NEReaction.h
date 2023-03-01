#pragma once

#include <xolotl/core/network/NETraits.h>
#include <xolotl/core/network/NucleationReaction.h>
#include <xolotl/core/network/ReSolutionReaction.h>

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
public:
	using Superclass =
		ProductionReaction<NEReactionNetwork, NEProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);
};

class NEDissociationReaction :
	public DissociationReaction<NEReactionNetwork, NEDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<NEReactionNetwork, NEDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
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

} // namespace network
} // namespace core
} // namespace xolotl
