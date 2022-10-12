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

private:
	// Alignment between loops.
	// Free and Trapped
	util::Array<double, 2> _align4by4{1.0, 0.3333333333333333};
	double _q4by4 = 0.0625;
	// Free, Trapped, and Junction
	util::Array<double, 3> _align4by7{
		1.0, 0.3333333333333333, 0.5773502691896258};
	double _q4by7 = 0.0357142857;
	// Free, Trapped, and Loop
	util::Array<double, 1> _align4by3{0.5773502691896258};
	double _q4by3 = 0.08333333333333333;
	// Junction and Junction
	util::Array<double, 4> _align7by7{
		1.0, 0.3333333333333333, 0.5773502691896258, 0.0};
	double _q7by7 = 0.02040816326530612;
	// Junction and Loop
	util::Array<double, 3> _align7by3{0.5773502691896258, 1.0, 0.0};
	double _q7by3 = 0.047619047619047616;
	// Loop and Loop
	util::Array<double, 2> _align3by3{1.0, 0.0};
	double _q3by3 = 0.1111111111111111;
};

class FeCrSinkReaction :
	public SinkReaction<FeCrReactionNetwork, FeCrSinkReaction>
{
public:
	using Superclass = SinkReaction<FeCrReactionNetwork, FeCrSinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();
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
