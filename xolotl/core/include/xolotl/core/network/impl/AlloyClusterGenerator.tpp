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
AlloyClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	result[0] = true;
	result[1] = true;
	result[2] = true;
	result[3] = true;
	result[4] = true;
	result[5] = true;

	int nAxis = (region[Species::V].begin() > 0) +
		(region[Species::I].begin() > 0) +
		(region[Species::PerfectV].begin() > 0) +
		(region[Species::FaultedV].begin() > 0) +
		(region[Species::FaultedI].begin() > 0) +
		(region[Species::PerfectI].begin() > 0);

	if (nAxis > 1) {
		result[0] = false;
		result[1] = false;
		result[2] = false;
		result[3] = false;
		result[4] = false;
		result[5] = false;
		return false;
	}

	if (nAxis == 0)
		return true;

	// I is always refined
	if (region[Species::I].begin() > 0)
		return true;

	// Smaller that the minimum size for grouping
	if (region[Species::V].begin() < _groupingMin &&
		region[Species::PerfectV].begin() < _groupingMin &&
		region[Species::FaultedV].begin() < _groupingMin &&
		region[Species::FaultedI].begin() < _groupingMin &&
		region[Species::PerfectI].begin() < _groupingMin) {
		return true;
	}

	// Too large
	if (region[Species::V].end() > _maxSize ||
		region[Species::PerfectV].end() > _maxSize ||
		region[Species::FaultedV].end() > _maxSize ||
		region[Species::FaultedI].end() > _maxSize ||
		region[Species::PerfectI].end() > _maxSize) {
		return true;
	}

	if (region[Species::V].begin() > 0 &&
		region[Species::V].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::V].begin(), 1) * 5.0e-2))
		result[0] = false;
	if (region[Species::PerfectV].begin() > 0 &&
		region[Species::PerfectV].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::PerfectV].begin(), 1) * 5.0e-2))
		result[1] = false;
	if (region[Species::FaultedV].begin() > 0 &&
		region[Species::FaultedV].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::FaultedV].begin(), 1) * 5.0e-2))
		result[2] = false;
	if (region[Species::FaultedI].begin() > 0 &&
		region[Species::FaultedI].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::FaultedI].begin(), 1) * 5.0e-2))
		result[5] = false;
	if (region[Species::PerfectI].begin() > 0 &&
		region[Species::PerfectI].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::PerfectI].begin(), 1) * 5.0e-2))
		result[4] = false;

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
AlloyClusterGenerator::select(const Region& region) const
{
	int nAxis = (region[Species::V].begin() > 0) +
		(region[Species::I].begin() > 0) +
		(region[Species::PerfectI].begin() > 0) +
		(region[Species::FaultedI].begin() > 0) +
		(region[Species::FaultedV].begin() > 0) +
		(region[Species::PerfectV].begin() > 0);

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
		if (region[Species::V].begin() > _maxSize)
			return false;

		// Perfect I
		if (region[Species::PerfectI].begin() > 0 &&
			region[Species::PerfectI].begin() <= _maxI)
			return false;
		if (region[Species::PerfectI].begin() > _maxSize)
			return false;

		// Faulted I
		if (region[Species::FaultedI].begin() > 0 &&
			region[Species::FaultedI].begin() <= _maxI)
			return false;
		if (region[Species::FaultedI].begin() > _maxSize)
			return false;

		// Faulted V
		if (region[Species::FaultedV].begin() > 0 &&
			region[Species::FaultedV].begin() <= _maxV)
			return false;
		if (region[Species::FaultedV].begin() > _maxSize)
			return false;

		// Perfect V
		if (region[Species::PerfectV].begin() > 0 &&
			region[Species::PerfectV].begin() <= _maxV)
			return false;
		if (region[Species::PerfectV].begin() > _maxSize)
			return false;
	}

	if (region[Species::PerfectV].begin() > 0 &&
		region[Species::PerfectV].end() - 1 <= _maxV)
		return false;

	if (region[Species::FaultedV].begin() > 0 &&
		region[Species::FaultedV].end() - 1 <= _maxV)
		return false;

	if (region[Species::PerfectI].begin() > 0 &&
		region[Species::PerfectI].end() - 1 <= _maxI)
		return false;

	if (region[Species::FaultedI].begin() > 0 &&
		region[Species::FaultedI].end() - 1 <= _maxI)
		return false;

	if (region[Species::V].begin() > _maxSize ||
		region[Species::PerfectV].begin() > _maxSize ||
		region[Species::FaultedV].begin() > _maxSize ||
		region[Species::PerfectI].begin() > _maxSize ||
		region[Species::FaultedI].begin() > _maxSize)
		return false;

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
AlloyClusterGenerator::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition lo(reg.getOrigin());
	double energy = 0.0;
	if (lo.isOnAxis(Species::PerfectI)) {
		for (auto j : makeIntervalRange(reg[Species::PerfectI])) {
			energy += 3.4 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::PerfectI].length();
	}
	if (lo.isOnAxis(Species::FaultedI)) {
		for (auto j : makeIntervalRange(reg[Species::FaultedI])) {
			energy += 3.4 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::FaultedI].length();
	}
	if (lo.isOnAxis(Species::FaultedV)) {
		for (auto j : makeIntervalRange(reg[Species::FaultedV])) {
			energy += 1.9 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::FaultedV].length();
	}
	if (lo.isOnAxis(Species::PerfectV)) {
		for (auto j : makeIntervalRange(reg[Species::PerfectV])) {
			energy += 1.9 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::PerfectV].length();
	}
	if (lo.isOnAxis(Species::V)) {
		for (auto j : makeIntervalRange(reg[Species::V])) {
			energy += 1.9 + 3.4 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::V].length();
	}
	if (lo.isOnAxis(Species::I)) {
		for (auto j : makeIntervalRange(reg[Species::I])) {
			energy += 3.4 + 3.5 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::I].length();
	}
	return 0.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
AlloyClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double migrationEnergy = util::infinity<double>;
	if (comp.isOnAxis(Species::V) and comp[Species::V] <= _maxV) {
		return 1.3;
	}
	if (comp.isOnAxis(Species::I)) {
		return 0.5;
	}
	return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
AlloyClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double diffusionFactor = 0.0;
	if (comp.isOnAxis(Species::V) and comp[Species::V] <= _maxV) {
		const double jumpDistance = latticeParameter / sqrt(2.0);
		constexpr double phononFrequency = 9.6e12;
		constexpr double jumpsPerPhonon = 1.0;
		constexpr double prefactorExponent = -1.0;
		return phononFrequency * jumpsPerPhonon * jumpDistance * jumpDistance *
			pow((double)comp[Species::V], prefactorExponent) / (6.0);
	}
	if (comp.isOnAxis(Species::I)) {
		const double jumpDistance = latticeParameter / sqrt(2.0);
		constexpr double phononFrequency = 9.6e12;
		constexpr double jumpsPerPhonon = 1.0;
		constexpr double prefactorExponent = -1.0;
		return phononFrequency * jumpsPerPhonon * jumpDistance * jumpDistance *
			pow((double)comp[Species::I], prefactorExponent) / (6.0);
	}
	return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
AlloyClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	const double prefactor =
		0.25 * latticeParameter * latticeParameter / ::xolotl::core::pi;
	const auto& reg = cluster.getRegion();
	Composition lo(reg.getOrigin());
	double radius = 0.0;
	if (lo.isOnAxis(Species::PerfectI)) {
		for (auto j : makeIntervalRange(reg[Species::PerfectI])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::perfectBurgers);
		}
		return radius / reg[Species::PerfectI].length();
	}
	if (lo.isOnAxis(Species::FaultedI)) {
		for (auto j : makeIntervalRange(reg[Species::FaultedI])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::faultedBurgers);
		}
		return radius / reg[Species::FaultedI].length();
	}
	if (lo.isOnAxis(Species::FaultedV)) {
		for (auto j : makeIntervalRange(reg[Species::FaultedV])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::faultedBurgers);
		}
		return radius / reg[Species::FaultedV].length();
	}
	if (lo.isOnAxis(Species::PerfectV)) {
		for (auto j : makeIntervalRange(reg[Species::PerfectV])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::perfectBurgers);
		}
		return radius / reg[Species::PerfectV].length();
	}
	if (lo.isOnAxis(Species::V)) {
		for (auto j : makeIntervalRange(reg[Species::V])) {
			radius += cbrt(0.75 * prefactor * latticeParameter * (double)j);
		}
		return radius / reg[Species::V].length();
	}
	if (lo.isOnAxis(Species::I)) {
		for (auto j : makeIntervalRange(reg[Species::I])) {
			radius += cbrt(0.75 * prefactor * latticeParameter * (double)j);
		}
		return radius / reg[Species::I].length();
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
