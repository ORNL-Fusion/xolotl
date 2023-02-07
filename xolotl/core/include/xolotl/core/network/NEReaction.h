#pragma once

#include <xolotl/core/network/FullReSolutionReaction.h>
#include <xolotl/core/network/NETraits.h>
#include <xolotl/core/network/NucleationReaction.h>
#include <xolotl/core/network/PartialReSolutionReaction.h>

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

class NEFullReSolutionReaction :
	public FullReSolutionReaction<NEReactionNetwork, NEFullReSolutionReaction>
{
public:
	using Superclass =
		FullReSolutionReaction<NEReactionNetwork, NEFullReSolutionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	void
	setSize();
};

class NEPartialReSolutionReaction :
	public PartialReSolutionReaction<NEReactionNetwork,
		NEPartialReSolutionReaction>
{
public:
	using Superclass = PartialReSolutionReaction<NEReactionNetwork,
		NEPartialReSolutionReaction>;

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
