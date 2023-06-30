#pragma once

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
	using IndexType = typename Superclass::IndexType;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);
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
	using IndexType = typename Superclass::IndexType;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
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
	using IndexType = typename Superclass::IndexType;

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
} // namespace network
} // namespace core
} // namespace xolotl
