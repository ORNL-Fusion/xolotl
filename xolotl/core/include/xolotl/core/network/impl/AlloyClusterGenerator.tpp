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
		(region[Species::Perfect].begin() > 0) +
		(region[Species::Frank].begin() > 0) +
		(region[Species::Faulted].begin() > 0) +
		(region[Species::Void].begin() > 0);

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

	// V, and I are always refined
	if (region[Species::V].begin() > 0)
		return true;
	if (region[Species::I].begin() > 0)
		return true;

	// Smaller than the minimum size for grouping
	if (region[Species::Void].begin() < _groupingMin &&
		region[Species::Faulted].begin() < _groupingMin &&
		region[Species::Frank].begin() < _groupingMin &&
		region[Species::Perfect].begin() < _groupingMin) {
		return true;
	}

	// Too large
	if (region[Species::Void].end() > _maxVoid && _maxVoid > 0)
		return true;
	if (region[Species::Faulted].end() > _maxSize && _maxSize > 0)
		return true;
	if (region[Species::Frank].end() > _maxSize && _maxSize > 0)
		return true;
	if (region[Species::Perfect].end() > _maxSize && _maxSize > 0)
		return true;

	if (region[Species::Void].begin() > 0 &&
		region[Species::Void].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::Void].begin(), 1) * 5.0e-2))
		result[1] = false;
	if (region[Species::Faulted].begin() > 0 &&
		region[Species::Faulted].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::Faulted].begin(), 1) * 5.0e-2))
		result[2] = false;
	if (region[Species::Frank].begin() > 0 &&
		region[Species::Frank].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::Frank].begin(), 1) * 5.0e-2))
		result[5] = false;
	if (region[Species::Perfect].begin() > 0 &&
		region[Species::Perfect].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::Perfect].begin(), 1) * 5.0e-2))
		result[4] = false;

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
AlloyClusterGenerator::select(const Region& region) const
{
	int nAxis = (region[Species::V].begin() > 0) +
		(region[Species::I].begin() > 0) +
		(region[Species::Perfect].begin() > 0) +
		(region[Species::Frank].begin() > 0) +
		(region[Species::Faulted].begin() > 0) +
		(region[Species::Void].begin() > 0);

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

		// Perfect
		if (region[Species::Perfect].begin() > 0 &&
			region[Species::Perfect].begin() < _maxI)
			return false;
		if (region[Species::Perfect].begin() > 0 &&
			region[Species::Perfect].begin() > _maxSize)
			return false;

		// Frank
		if (region[Species::Frank].begin() > 0 &&
			region[Species::Frank].begin() <= _maxI)
			return false;
		if (region[Species::Frank].begin() > 0 &&
			region[Species::Frank].begin() > _maxSize)
			return false;

		// Faulted
		if (region[Species::Faulted].begin() > 0 &&
			region[Species::Faulted].begin() <= _maxV)
			return false;
		if (region[Species::Faulted].begin() > 0 &&
			region[Species::Faulted].begin() > _maxSize)
			return false;

		// Void
		if (region[Species::Void].begin() > 0 &&
			region[Species::Void].begin() <= _maxV)
			return false;
		if (region[Species::Void].begin() > 0 &&
			region[Species::Void].begin() > _maxVoid)
			return false;
	}

	if (region[Species::V].begin() == 0 && region[Species::I].begin() == 0 &&
		region[Species::Void].end() - 1 <= _maxV &&
		region[Species::Faulted].end() - 1 <= _maxV &&
		region[Species::Perfect].end() - 1 <= _maxI &&
		region[Species::Frank].end() - 1 <= _maxI)
		return false;

	if (region[Species::Void].begin() > _maxVoid ||
		region[Species::Faulted].begin() > _maxSize ||
		region[Species::Perfect].begin() > _maxSize ||
		region[Species::Frank].begin() > _maxSize)
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
	if (lo.isOnAxis(Species::Perfect)) {
		for (auto j : makeIntervalRange(reg[Species::Perfect])) {
			energy += 3.4 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::Perfect].length();
	}
	if (lo.isOnAxis(Species::Frank)) {
		for (auto j : makeIntervalRange(reg[Species::Frank])) {
			energy += 3.4 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::Frank].length();
	}
	if (lo.isOnAxis(Species::Faulted)) {
		for (auto j : makeIntervalRange(reg[Species::Faulted])) {
			energy += 1.9 + 2.0 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::Faulted].length();
	}
	if (lo.isOnAxis(Species::Void)) {
		for (auto j : makeIntervalRange(reg[Species::Void])) {
			energy += 1.9 + 3.4 * (pow((double)j, 2.0 / 3.0) - 1.0);
		}
		return energy / reg[Species::Void].length();
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
	if (comp.isOnAxis(Species::Perfect)) {
		return 0.5;
	}
	if (comp.isOnAxis(Species::V)) {
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
	if (comp.isOnAxis(Species::Perfect)) {
		if (comp[Species::Perfect] < 70) {
			const double jumpDistance = latticeParameter / sqrt(2.0);
			constexpr double phononFrequency = 9.6e12;
			constexpr double jumpsPerPhonon = 1.0;
			constexpr double prefactorExponent = -1.0;
			return phononFrequency * jumpsPerPhonon * jumpDistance *
				jumpDistance *
				pow((double)comp[Species::Perfect], prefactorExponent) / (6.0);
		}
		else
			return diffusionFactor;
	}
	if (comp.isOnAxis(Species::V)) {
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
	if (lo.isOnAxis(Species::Perfect)) {
		for (auto j : makeIntervalRange(reg[Species::Perfect])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::perfectBurgers);
		}
		return radius / reg[Species::Perfect].length();
	}
	if (lo.isOnAxis(Species::Frank)) {
		for (auto j : makeIntervalRange(reg[Species::Frank])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::frankBurgers);
		}
		return radius / reg[Species::Frank].length();
	}
	if (lo.isOnAxis(Species::Faulted)) {
		for (auto j : makeIntervalRange(reg[Species::Faulted])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::faultedBurgers);
		}
		return radius / reg[Species::Faulted].length();
	}
	if (lo.isOnAxis(Species::Void)) {
		for (auto j : makeIntervalRange(reg[Species::Void])) {
			radius += cbrt(0.75 * prefactor * latticeParameter * (double)j);
		}
		return radius / reg[Species::Void].length();
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
