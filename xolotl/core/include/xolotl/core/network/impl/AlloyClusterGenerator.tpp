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
	// Declarations
	using alloy::getMaxHePerV;
	using detail::toIndex;

	// Initialize the result array that indicates which axis needs to be refined
	// more
	for (auto& r : result) {
		r = true;
	}

	// Check if others begin at 0
	auto othersBeginAtZero = [](const Region& reg, Species species) {
		for (auto s : NetworkType::getSpeciesRange()) {
			if (s.value != species && reg[s].begin() != 0) {
				return false;
			}
		}
		return true;
	};

	// Count the number of axis it is on
	int nAxis = 0;
	for (auto s : NetworkType::getSpeciesRange()) {
		if (region[s].begin() > 0) {
			nAxis++;
		}
	}

	// Cannot be on more than 2 axis
	if (nAxis > 2) {
		for (auto& r : result) {
			r = false;
		}
		return false;
	}

	// Cannot be 0 axis
	if (nAxis == 0)
		return true;

	Composition lo = region.getOrigin();
	Composition hi = region.getUpperLimitPoint();

	// Smaller that the minimum size for grouping
	if (region[Species::He].begin() < _groupingMin &&
		region[Species::V].begin() < _groupingMin &&
		region[Species::PerfectV].begin() < _groupingMin &&
		region[Species::FaultedV].begin() < _groupingMin &&
		region[Species::FaultedI].begin() < _groupingMin &&
		region[Species::PerfectI].begin() < _groupingMin) {
		return true;
	}

	// I is always refined
	if (region[Species::I].begin() > 0)
		return true;

	// Grouping loops
	if (lo[Species::PerfectV] > 0) {
		if (lo[Species::PerfectV] < _groupingMin &&
			othersBeginAtZero(region, Species::PerfectV)) {
			return true;
		}
		if (region[Species::PerfectV].end() > _maxLoopSize) {
			return true;
		}
		if (region[Species::PerfectV].length() <
			util::max((double)(_groupingWidthB + 1),
				pow(region[Species::PerfectV].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::PerfectV)] = false;
		}
	}

	if (lo[Species::FaultedV] > 0) {
		if (lo[Species::FaultedV] < _groupingMin &&
			othersBeginAtZero(region, Species::FaultedV)) {
			return true;
		}
		if (region[Species::FaultedV].end() > _maxLoopSize) {
			return true;
		}
		if (region[Species::FaultedV].length() <
			util::max((double)(_groupingWidthB + 1),
				pow(region[Species::FaultedV].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::FaultedV)] = false;
		}
	}

	if (lo[Species::PerfectI] > 0) {
		if (lo[Species::PerfectI] < _groupingMin &&
			othersBeginAtZero(region, Species::PerfectI)) {
			return true;
		}
		if (region[Species::PerfectI].end() > _maxLoopSize) {
			return true;
		}
		if (region[Species::PerfectI].length() <
			util::max((double)(_groupingWidthB + 1),
				pow(region[Species::PerfectI].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::PerfectI)] = false;
		}
	}

	if (lo[Species::FaultedI] > 0) {
		if (lo[Species::FaultedI] < _groupingMin &&
			othersBeginAtZero(region, Species::FaultedI)) {
			return true;
		}
		if (region[Species::FaultedI].end() > _maxLoopSize) {
			return true;
		}
		if (region[Species::FaultedI].length() <
			util::max((double)(_groupingWidthB + 1),
				pow(region[Species::FaultedI].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::FaultedI)] = false;
		}
	}

	// Group pure V
	if (lo[Species::V] > 0 and othersBeginAtZero(region, Species::V)) {
		if (region[Species::V].end() > _maxSize) {
			result[toIndex(Species::V)] = true;
			return true;
		}
		if (region[Species::V].length() <
			util::max((double)(_groupingWidthB + 1),
				pow(region[Species::V].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::V)] = false;
		}
	}

	// He is never grouped
	if (hi[Species::He] > 1 && othersBeginAtZero(region, Species::He)) {
		return true;
	}

	// Else refine around the edge
	if (hi[Species::V] > 1) {
		if (lo[Species::He] <= getMaxHePerV(hi[Species::V] - 1, _hevRatio) &&
			hi[Species::He] - 1 >=
				getMaxHePerV(lo[Species::V] - 1, _hevRatio)) {
			return true;
		}
	}

	// HeV phasespace
	if (region[Species::V].length() <
		util::max((double)(_groupingWidthA + 1),
			pow(region[Species::V].begin(), 0.75) * 0.5)) {
		result[toIndex(Species::V)] = false;
	}

	if (region[Species::He].length() <
		util::max((double)(_groupingWidthA + 1),
			pow(region[Species::He].begin(), 0.75) * 0.5)) {
		result[toIndex(Species::He)] = false;
	}

	// Border case for He and V
	if (hi[Species::V] > _maxSize) {
		result[toIndex(Species::V)] = true;
	}
	if (lo[Species::V] == 0) {
		result[toIndex(Species::V)] = true;
	}

	if (hi[Species::He] > _maxSize) {
		result[toIndex(Species::He)] = true;
	}
	if (lo[Species::He] == 0) {
		result[toIndex(Species::He)] = true;
	}

	int axis = 0;
	for (auto& r : result) {
		axis += r;
	}

	if (axis == 0)
		return false;

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
AlloyClusterGenerator::select(const Region& region) const
{
	// Declarations
	using alloy::getMaxHePerV;
	using detail::toIndex;

	// Remove 0
	auto isZeroPoint = [](const Region& reg) {
		for (const auto& ival : reg) {
			if (ival.end() != 1) {
				return false;
			}
		}
		return true;
	};
	if (isZeroPoint(region)) {
		return false;
	}

	auto othersEndAtOne = [](const Region& reg, Species species) {
		for (auto s : NetworkType::getSpeciesRange()) {
			if (s.value != species && reg[s].end() != 1) {
				return false;
			}
		}
		return true;
	};

	Composition lo = region.getOrigin();
	Composition hi = region.getUpperLimitPoint();

	// Interstitials
	if (region[Species::I].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::I)) {
		return false;
	}
	if (region[Species::I].begin() > _maxI) {
		return false;
	}

	// Helium
	if (region[Species::He].begin() > 1 &&
		othersEndAtOne(region, Species::He)) {
		return false;
	}
	// He can only cluster with V
	if (region[Species::He].begin() > 1 && hi[Species::V] == 1) {
		return false;
	}

	// Loops
	if (region[Species::PerfectV].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::PerfectV)) {
		return false;
	}
	if (region[Species::PerfectV].begin() > 0 && region.isSimplex() &&
		region[Species::PerfectV].begin() <= _maxV) {
		return false;
	}
	if (region[Species::PerfectV].begin() > _maxLoopSize) {
		return false;
	}

	if (region[Species::FaultedV].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::FaultedV)) {
		return false;
	}
	if (region[Species::FaultedV].begin() > 0 && region.isSimplex() &&
		region[Species::FaultedV].begin() <= _maxV) {
		return false;
	}
	if (region[Species::FaultedV].begin() > _maxLoopSize) {
		return false;
	}

	if (region[Species::PerfectI].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::PerfectI)) {
		return false;
	}
	if (region[Species::PerfectI].begin() > 0 && region.isSimplex() &&
		region[Species::PerfectI].begin() <= _maxI) {
		return false;
	}
	if (region[Species::PerfectI].begin() > _maxLoopSize) {
		return false;
	}

	if (region[Species::FaultedI].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::FaultedI)) {
		return false;
	}
	if (region[Species::FaultedI].begin() > 0 && region.isSimplex() &&
		region[Species::FaultedI].begin() <= _maxI) {
		return false;
	}
	if (region[Species::FaultedI].begin() > _maxLoopSize) {
		return false;
	}
	if (region[Species::PerfectV].end() - 1 <= _maxV &&
		region[Species::FaultedV].end() - 1 <= _maxV &&
		region[Species::PerfectI].end() - 1 <= _maxI &&
		region[Species::FaultedI].end() - 1 <= _maxI &&
		region[Species::He].begin() == 0 && region[Species::V].begin() == 0 &&
		region[Species::I].begin() == 0) {
		return false;
	}

	// Vacancy
	if (region[Species::V].begin() > _maxSize) {
		return false;
	}

	// The edge
	if (region[Species::V].end() > 1) {
		auto hiV = util::min(hi[Species::V] - 1, _maxSize);
		auto hiHe =
			util::min(hi[Species::He] - 1, getMaxHePerV(_maxSize, _hevRatio));

		// Too many helium
		if (lo[Species::He] >
			util::min(getMaxHePerV(hiV, _hevRatio), _maxSize)) {
			return false;
		}
	}

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
AlloyClusterGenerator::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// These are actually never called
	const auto& reg = cluster.getRegion();
	Composition lo(reg.getOrigin());
	double energy = 0.0;
	if (lo[Species::He] == 1) {
		for (auto j : makeIntervalRange(reg[Species::He])) {
			for (auto i : makeIntervalRange(reg[Species::V])) {
				double ratio = ((double)j) / (double)i;
				energy += 3.2 * (1.0 - exp(-0.4167 * pow(ratio, 0.9477)));
			}
		}
		return energy / (reg[Species::He].length() * reg[Species::V].length());
	}
	if (lo[Species::He] > 1) {
		for (auto j : makeIntervalRange(reg[Species::He])) {
			for (auto i : makeIntervalRange(reg[Species::V])) {
				double ratio_before = ((double)j - 1.0) / (double)i;
				double ratio = ((double)j) / (double)i;
				energy += (2.101 + 1.545 * (ratio_before / ratio)) * 3.2 *
					(1.0 - exp(-0.4167 * pow(ratio, 0.9477)));
			}
		}
		return energy / (reg[Species::He].length() * reg[Species::V].length());
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
	if (comp.isOnAxis(Species::He)) {
		return 0.13;
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
	if (comp.isOnAxis(Species::He)) {
		return 1.0e11;
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
	if (lo[Species::V] > 0) {
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
	if (lo.isOnAxis(Species::He)) {
		return impurityRadius;
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
