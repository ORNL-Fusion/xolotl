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
TransformReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex, double)
{
	auto cl = this->_clusterData->getCluster(_reactant);

	auto exponent = this->asDerived()->getExponent();
	auto barrier = this->asDerived()->getBarrier();
	auto size = this->asDerived()->getSize();

	double strength = 1.0e13 *
		exp(-barrier /
			(::xolotl::core::kBoltzmann *
				this->_clusterData->temperature(gridIndex))) /
		pow(size + 1.0, exponent);

	return strength;
}
} // namespace network
} // namespace core
} // namespace xolotl
