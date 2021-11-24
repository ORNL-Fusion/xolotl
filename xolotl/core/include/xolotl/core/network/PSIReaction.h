#pragma once

#include <xolotl/core/network/BurstingReaction.h>
#include <xolotl/core/network/PSITraits.h>
#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/TrapMutationReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIProductionReaction :
	public ProductionReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIProductionReaction<TSpeciesEnum>>
{
public:
	using Superclass = ProductionReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIProductionReaction<TSpeciesEnum>>;

	using Superclass::Superclass;
};

template <typename TSpeciesEnum>
class PSIDissociationReaction :
	public DissociationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIDissociationReaction<TSpeciesEnum>>
{
public:
	using Superclass = DissociationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIDissociationReaction<TSpeciesEnum>>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy();
};

template <typename TSpeciesEnum>
class PSISinkReaction :
	public SinkReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSISinkReaction<TSpeciesEnum>>
{
public:
	using Superclass = SinkReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSISinkReaction<TSpeciesEnum>>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();
};

template <typename TSpeciesEnum>
class PSITrapMutationReaction :
	public TrapMutationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSITrapMutationReaction<TSpeciesEnum>>
{
public:
	using Superclass = TrapMutationReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSITrapMutationReaction<TSpeciesEnum>>;
	using Superclass::Superclass;
};
template <typename TSpeciesEnum>
class PSIBurstingReaction :
	public BurstingReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIBurstingReaction<TSpeciesEnum>>
{
public:
	using Superclass = BurstingReaction<PSIReactionNetwork<TSpeciesEnum>,
		PSIBurstingReaction<TSpeciesEnum>>;
	using IndexType = typename Superclass::IndexType;
	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getAppliedRate(IndexType gridIndex) const;
};
} // namespace network
} // namespace core
} // namespace xolotl
