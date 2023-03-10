#pragma once

#include <xolotl/core/network/FeCrTraits.h>
#include <xolotl/core/network/SinkReaction.h>
#include <xolotl/core/network/TransformReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class FeCrReactionNetwork;

class FeCrProductionReaction :
	public ProductionReaction<FeCrReactionNetwork, FeCrProductionReaction>
{
public:
	using Superclass =
		ProductionReaction<FeCrReactionNetwork, FeCrProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeNetSigma(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex);
};

class FeCrDissociationReaction :
	public DissociationReaction<FeCrReactionNetwork, FeCrDissociationReaction>
{
public:
	using Superclass =
		DissociationReaction<FeCrReactionNetwork, FeCrDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy(double time = 0.0);
};

class FeCrSinkReaction :
	public SinkReaction<FeCrReactionNetwork, FeCrSinkReaction>
{
public:
	using Superclass = SinkReaction<FeCrReactionNetwork, FeCrSinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();

	KOKKOS_INLINE_FUNCTION
	double
	computeNetSigma(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex);
};

class FeCrTransformReaction :
	public TransformReaction<FeCrReactionNetwork, FeCrTransformReaction>
{
public:
	using Superclass =
		TransformReaction<FeCrReactionNetwork, FeCrTransformReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSize();

	KOKKOS_INLINE_FUNCTION
	double
	getExponent();

	KOKKOS_INLINE_FUNCTION
	double
	getBarrier();
};
} // namespace network
} // namespace core
} // namespace xolotl
