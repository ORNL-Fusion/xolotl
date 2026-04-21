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
	result[1] = true;

	// H is never grouped
	if (region[Species::H].begin() > 0) {
		return true;
	}

	// V is never grouped
	if (region[Species::V].begin() > 0) {
		return true;
	}

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
LiClusterGenerator::select(const Region& region) const
{
        int nAxis = (region[Species::H].begin() > 0) +
		(region[Species::V].begin() > 0);

	if (nAxis > 1) {
		return false;
	}

	if (region.isSimplex()) {
		// Each cluster should be on one axis and one axis only
		if (nAxis != 1) {
			return false;
		}
	}

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
LiClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double migrationEnergy = util::infinity<double>;
	if (comp.isOnAxis(Species::H)) {
		return 0.5;
	}
	
	return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
LiClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double diffusionFactor = 0.0;
	if (comp.isOnAxis(Species::H)) {
	        return 4.07e15;
	}
	
	return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
LiClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	if (comp.isOnAxis(Species::H)) {
	        return 0.0529;
	}
	
	return latticeParameter;
}
} // namespace network
} // namespace core
} // namespace xolotl
