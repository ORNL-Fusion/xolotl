#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
KOKKOS_INLINE_FUNCTION
bool
LiClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	result[0] = true;

	// H is never grouped
	if (region[Species::H].begin() > 0) {
		return true;
	}

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
LiClusterGenerator::select(const Region& region) const
{
	// Remove 0
	if (region[Species::H].end() == 1) {
		return false;
	}

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
LiClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// From paper: 141.5 kJ mol-1
	//	return 0.0;
	return 141500.0 * xolotl::core::kBoltzmann / 8.314;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
LiClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	// From paper: 1.29e-5 m2 s-1
	return 1.29e13;
	//	return 6.86e6;
	//	return 1.18e6;
	//	return 1.0e4;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
LiClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	// TODO
	return impurityRadius;
}
} // namespace network
} // namespace core
} // namespace xolotl
