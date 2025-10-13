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
RPVClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	for (auto& r : result) {
		r = true;
	}

	// I is never grouped
	if (region[Species::I].begin() > 0) {
		return true;
	}

	// V is grouped on its own
	if (region[Species::V].end() > 1 && region[Species::I].begin() == 0) {
		if (region[Species::V].begin() < _groupingMin)
			return true;
		if (region[Species::V].end() > _maxV) {
			return true;
		}
		if (region[Species::V].begin() > 0 &&
			region[Species::V].length() <
				util::max((double)(_groupingWidthV + 1),
					pow(region[Species::V].begin(), 0.75) * 0.1))
			result[0] = false;
		else
			return true;
	}

	// Edges
	if (region[Species::V].end() > _maxV + 1) {
		return true;
	}

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
RPVClusterGenerator::select(const Region& region) const
{
	// Remove 0
	if (region[Species::V].end() == 1 && region[Species::I].end() == 1) {
		return false;
	}

	// Interstitials
	if (region[Species::I].begin() > 0 && region[Species::V].begin() > 0) {
		return false;
	}

	// Vacancy
	if (region[Species::V].begin() > _maxV && region[Species::I].end() == 1) {
		return false;
	}
	if (region[Species::V].begin() > _maxV) {
		return false;
	}

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
RPVClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// I migration energy in eV
	constexpr double iOneMigrationEnergy = 0.22;
	// V migration energies in eV
	constexpr Kokkos::Array<double, 5> vMigration = {
		0.0, 0.63, 0.50, 0.36, 0.44};

	const auto& reg = cluster.getRegion();
	double migrationEnergy = util::infinity<double>;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			if (comp[Species::I] == 1) {
				migrationEnergy = iOneMigrationEnergy;
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto amtV = comp[Species::V];
			if (amtV < vMigration.size()) {
				migrationEnergy = vMigration[amtV];
			}
		}
	}
	return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
RPVClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	// I diffusion factors in nm^2/s
	constexpr double iOneDiffusionFactor = 2.8e+10;
	// V diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 5> vDiffusion = {
		0.0, 8.2e+11, 4.1e+11, 2.73e+11, 2.05e+11};

	const auto& reg = cluster.getRegion();
	double diffusionFactor = 0.0;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			if (comp[Species::I] == 1) {
				diffusionFactor = iOneDiffusionFactor;
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto amtV = comp[Species::V];
			if (amtV < vDiffusion.size()) {
				diffusionFactor = vDiffusion[amtV];
			}
		}
	}

	return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
RPVClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	const auto& reg = cluster.getRegion();
	double radius = 0.0;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			radius = latticeParameter * cbrt(3.0 / ::xolotl::core::pi) * 0.5;
		}
		else {
			radius = latticeParameter *
				cbrt((3.0 * comp[Species::V]) / ::xolotl::core::pi) * 0.5;
		}
	}
	else {
		// Loop on the V range
		for (auto j : makeIntervalRange(reg[Species::V])) {
			radius += latticeParameter *
				cbrt((3.0 * (double)j) / ::xolotl::core::pi) * 0.5;
		}
		// Average the radius
		radius /= reg[Species::V].length();
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
