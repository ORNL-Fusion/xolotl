#pragma once

#include <plsm/EnumIndexed.h>

#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
SinkReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex, double time)
{
	auto cl = this->_clusterData->getCluster(_reactant);
	double dc = cl.getDiffusionCoefficient(gridIndex);

	double strength = this->asDerived()->getSinkBias() *
		this->asDerived()->getSinkStrength() * dc;

	return strength;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
SinkReaction<TNetwork, TDerived>::computeNetSigma(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// Check if our cluster is on the left side of this reaction
	if (clusterId == _reactant) {
		return this->asDerived()->getSinkBias() *
			this->asDerived()->getSinkStrength();
	}

	// This cluster is not part of the reaction
	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
