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
ZrClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	for (auto& r : result) {
		r = true;
	}

	int nAxis = (region[Species::V].begin() > 0) +
		(region[Species::I].begin() > 0) + (region[Species::Basal].begin() > 0);

	if (nAxis > 1) {
		for (auto& r : result) {
			r = false;
		}
		return false;
	}

	// Smaller that the minimum size for grouping
	if (region[Species::V].begin() < _groupingMin &&
		region[Species::Basal].begin() < _groupingMin &&
		region[Species::I].begin() < _groupingMin) {
		return true;
	}

	// Too large
	if (region[Species::V].end() > _maxV ||
		region[Species::Basal].end() > _maxB ||
		region[Species::I].end() > _maxI) {
		return true;
	}

	auto loV = region[Species::V].begin();
	if (loV > 0 &&
		region[Species::V].length() <
			util::max((double)(_groupingWidth + 1), loV * 2.0e-2))
		result[0] = false;
	auto loB = region[Species::Basal].begin();
	if (loB > 0 &&
		region[Species::Basal].length() <
			util::max((double)(_groupingWidth + 1), loB * 2.0e-2))
		result[1] = false;
	auto loI = region[Species::I].begin();
	if (loI > 0 &&
		region[Species::I].length() <
			util::max((double)(_groupingWidth + 1), loI * 2.0e-2))
		result[2] = false;

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
ZrClusterGenerator::select(const Region& region) const
{
	int nAxis = (region[Species::V].begin() > 0) +
		(region[Species::I].begin() > 0) + (region[Species::Basal].begin() > 0);

	if (nAxis > 1) {
		return false;
	}

	if (region.isSimplex()) {
		// Each cluster should be on one axis and one axis only
		if (nAxis != 1) {
			return false;
		}

		// I
		if (region[Species::I].begin() > _maxI)
			return false;

		// V
		if (region[Species::V].begin() > _maxV)
			return false;

		// Basal
		if (region[Species::Basal].begin() > _maxB)
			return false;
	}

	if (region[Species::V].begin() > _maxV ||
		region[Species::Basal].begin() > _maxB ||
		region[Species::I].begin() > _maxI)
		return false;

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
ZrClusterGenerator::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	return 0.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
ZrClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	constexpr double defaultMigration = -1.0;
	constexpr double iNineMigration = 0.10;

	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double migrationEnergy = util::infinity<double>;

	if (comp.isOnAxis(Species::V) && comp[Species::V] <= 6) {
		return defaultMigration;
	}
	if (comp.isOnAxis(Species::I)) {
		if (comp[Species::I] <= 3)
			return defaultMigration;
		if (comp[Species::I] == 9)
			return iNineMigration;
	}

	return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
ZrClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	constexpr double defaultDiffusion = 1.0;
	// constexpr double iNineDiffusion = 2.2e+11;
	constexpr double iNineDiffusion = 0.0;

	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double diffusionFactor = 0.0;

	if (comp.isOnAxis(Species::V) && comp[Species::V] <= 6) {
		return defaultDiffusion;
	}

	if (comp.isOnAxis(Species::I)) {
		if (comp[Species::I] <= 3)
			return defaultDiffusion;
		if (comp[Species::I] == 9)
			return iNineDiffusion;
	}

	return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
ZrClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition lo(reg.getOrigin());
	double radius = 0.0;
	int basalTransitionSize = 91;

	// jmr: rn = (3nOmega/4pi)^1/3 [nm] for n < 10
	// jmr: Note that (3Omega/4pi) = 5.586e-3 nm^3, where Omega = 0.0234 nm^3
	// jmr: For prismatic loops (n > 9): rn = 0.163076*sqrt(n) [nm]
	// jmr: For basal loops (n > 9): rn = 0.169587*sqrt(n) [nm]

	if (lo.isOnAxis(Species::V)) {
		for (auto j : makeIntervalRange(reg[Species::V])) {
			if (lo[Species::V] < 10)
				radius += pow(5.586e-3 * (double)j, 1.0 / 3.0);
			else
				radius += 0.163076 * pow((double)j, 0.5);
		}
		return radius / reg[Species::V].length();
	}

	// adding basal
	if (lo.isOnAxis(Species::Basal)) {
		for (auto j : makeIntervalRange(reg[Species::Basal])) {
			// Treat the case for faulted basal pyramids
			// Estimate a spherical radius based on equivalent surface area
			if (lo[Species::Basal] < basalTransitionSize) {
				double Sb = sqrt(3.0) / 2.0 * 3.232 * 3.232 *
					(double)j; // Basal surface area
				double Sp = 3.232 / 2.0 *
					sqrt(3.0 * 3.232 * 3.232 + 4.0 * 5.17 * 5.17) *
					(double)j; // Prismatic surface area
				radius += sqrt((Sb + Sp) / (4.0 * pi)) / 10.0;
			}

			// Treat the case of a basal c-loop
			else
				radius += 0.169587 * sqrt((double)j);
		}
		return radius / reg[Species::Basal].length();
	}

	if (lo.isOnAxis(Species::I)) {
		for (auto j : makeIntervalRange(reg[Species::I])) {
			if (lo[Species::I] < 10)
				radius += pow(5.586e-3 * (double)j, 1.0 / 3.0);
			else
				radius += 0.163076 * pow((double)j, 0.5);
		}
		return radius / reg[Species::I].length();
	}
	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
