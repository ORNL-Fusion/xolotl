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
FeCrClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	for (auto& r : result) {
		r = true;
	}

	int nAxis = 0;
	for (auto s : NetworkType::getSpeciesRange()) {
		if (region[s].begin() > 0) {
			nAxis++;
		}
	}

	if (nAxis > 2) {
		for (auto& r : result) {
			r = false;
		}
		return false;
	}

	//	if (nAxis == 2) {
	//		if (region[Species::Trap].end() - 1 < 1) {
	//			for (auto& r : result) {
	//				r = false;
	//			}
	//			return false;
	//		}
	//	}
	//
	//	if (nAxis == 0)
	//		return true;
	//
	//	// I and Complex are always refined
	//	if (region[Species::I].begin() > 0)
	//		return true;
	//	if (region[Species::Complex].begin() > 0)
	//		return true;
	//
	//	// Smaller that the minimum size for grouping
	//	if (region[Species::V].begin() < _groupingMin ||
	//		region[Species::Free].begin() < _groupingMin ||
	//		region[Species::Loop].begin() < _groupingMin ||
	//		region[Species::Trapped].begin() < _groupingMin ||
	//		region[Species::Junction].begin() < _groupingMin) {
	//		return true;
	//	}

	// Too large
	if (region[Species::V].begin() > _maxSize &&
		region[Species::Free].begin() > _maxSize &&
		region[Species::Trapped].begin() > _maxSize &&
		region[Species::Junction].begin() > _maxSize &&
		region[Species::Loop].begin() > _maxSize) {
		for (auto& r : result) {
			r = false;
		}
		return false;
	}

	//	if (region[Species::V].end() > _maxSize ||
	//		region[Species::Free].end() > _maxSize ||
	//		region[Species::Trapped].end() > _maxSize ||
	//		region[Species::Junction].end() > _maxSize ||
	//		region[Species::Loop].end() > _maxSize) {
	//		return true;
	//	}

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
FeCrClusterGenerator::select(const Region& region) const
{
	int nAxis = 0;
	for (auto s : NetworkType::getSpeciesRange()) {
		if (region[s].begin() > 0) {
			nAxis++;
		}
	}

	// V, I, Free are on one axis
	if (region[Species::I].begin() > 0 || region[Species::V].begin() > 0 ||
		region[Species::Free].begin() > 0) {
		if (nAxis > 1) {
			return false;
		}

		// I
		if (region[Species::I].begin() > _maxI)
			return false;

		// V
		if (region[Species::V].begin() > _maxSize)
			return false;

		// Free
		if (region[Species::Free].begin() > 0 &&
			region[Species::Free].begin() <= _maxI)
			return false;
		if (region[Species::Free].begin() > _maxSize)
			return false;
	}
	else {
		// Trapped, Junction, Complex, Loop are also on Trap
		if (nAxis > 2) {
			return false;
		}

		// Each region should have Trap = 1
		if (region[Species::Trapped].begin() > 0 ||
			region[Species::Junction].begin() > 0 ||
			region[Species::Complex].begin() > 0 ||
			region[Species::Loop].begin() > 0) {
			if (region[Species::Trap].begin() != 1)
				return false;
		}

		// Complex
		if (region[Species::Complex].begin() > _maxI)
			return false;

		// Junction
		if (region[Species::Junction].begin() > 0 &&
			region[Species::Junction].end() - 1 < _minJunction)
			return false;
		if (region[Species::Junction].begin() > _maxSize)
			return false;

		// Trapped
		if (region[Species::Trapped].begin() > 0 &&
			region[Species::Trapped].end() - 1 <= _maxI)
			return false;
		if (region[Species::Trapped].begin() > _maxSize)
			return false;

		// Loop
		if (region[Species::Loop].begin() > 0 &&
			region[Species::Loop].end() - 1 <= _maxI)
			return false;
		if (region[Species::Loop].begin() > _maxSize)
			return false;
	}

	if (region.isSimplex()) {
		// Remove 0
		if (nAxis == 0)
			return false;
	}

	//	if (region[Species::V].begin() > _maxSize ||
	//		region[Species::Free].begin() > _maxSize ||
	//		region[Species::Trapped].begin() > _maxSize ||
	//		region[Species::Complex].begin() > _maxSize ||
	//		region[Species::Junction].begin() > _maxSize ||
	//		region[Species::Loop].begin() > _maxSize)
	//		return false;

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
FeCrClusterGenerator::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition lo(reg.getOrigin());
	double energy = 0.0;
	if (lo.isOnAxis(Species::V) && lo[Species::V] == 1) {
		return 1.73;
	}
	return 0.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
FeCrClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double migrationEnergy = util::infinity<double>;
	if (comp.isOnAxis(Species::Free)) {
		return 0.1;
	}
	if (comp.isOnAxis(Species::V)) {
		switch (comp[Species::V]) {
		case 1:
			return 0.67;
		case 2:
			return 0.62;
		case 3:
			return 0.35;
		case 4:
			return 0.48;
		default:
			return migrationEnergy;
		}
	}
	if (comp.isOnAxis(Species::I)) {
		switch (comp[Species::I]) {
		case 1:
			return 0.34;
		case 2:
			return 0.42;
		case 3:
			return 0.43;
		default:
			return migrationEnergy;
		}
	}
	return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
FeCrClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double diffusionFactor = 0.0;
	// FeCr phonon frequency
	constexpr double phononFrequency = ::xolotl::core::fecrPhononFrequency;
	if (comp.isOnAxis(Species::Free)) {
		const double jumpDistance = latticeParameter * sqrt(3.0) / 2.0;
		constexpr double prefactorExponent = -0.7;
		return phononFrequency * jumpDistance * jumpDistance *
			pow((double)comp[Species::Free], prefactorExponent) / (2.0);
	}
	if (comp.isOnAxis(Species::V)) {
		if (comp[Species::V] < 5) {
			const double jumpDistance = latticeParameter * sqrt(3.0) / 2.0;
			constexpr double prefactorExponent = -1.0;
			return phononFrequency * jumpDistance * jumpDistance *
				pow((double)comp[Species::V], prefactorExponent) / (6.0);
		}
	}
	if (comp.isOnAxis(Species::I)) {
		const double jumpDistance = latticeParameter * sqrt(3.0) / 2.0;
		constexpr double prefactorExponent = -1.0;
		return phononFrequency * jumpDistance * jumpDistance *
			pow((double)comp[Species::I], prefactorExponent) / (6.0);
	}
	return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
FeCrClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	const double prefactor =
		0.5 * latticeParameter * latticeParameter / ::xolotl::core::pi;
	const auto& reg = cluster.getRegion();
	Composition lo(reg.getOrigin());
	double radius = 0.0;
	if (lo.isOnAxis(Species::Free)) {
		for (auto j : makeIntervalRange(reg[Species::Free])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::fecrBurgers);
		}
		return radius / reg[Species::Free].length();
	}
	if (lo[Species::Trapped] > 0) {
		for (auto j : makeIntervalRange(reg[Species::Trapped])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::fecrBurgers);
		}
		return radius / reg[Species::Trapped].length();
	}
	if (lo[Species::Junction] > 0) {
		for (auto j : makeIntervalRange(reg[Species::Junction])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::fecrBurgers);
		}
		return radius / reg[Species::Junction].length();
	}
	if (lo[Species::Loop] > 0) {
		for (auto j : makeIntervalRange(reg[Species::Loop])) {
			radius +=
				sqrt(((double)j * prefactor) / ::xolotl::core::fecrLoopBurgers);
		}
		return radius / reg[Species::Loop].length();
	}
	if (lo.isOnAxis(Species::V)) {
		for (auto j : makeIntervalRange(reg[Species::V])) {
			radius +=
				pow(0.75 * prefactor * latticeParameter * (double)j, 1.0 / 3.0);
		}
		return radius / reg[Species::V].length();
	}
	if (lo.isOnAxis(Species::I)) {
		for (auto j : makeIntervalRange(reg[Species::I])) {
			radius +=
				pow(0.75 * prefactor * latticeParameter * (double)j, 1.0 / 3.0);
			;
		}
		return radius / reg[Species::I].length();
	}
	if (lo[Species::Complex] > 0) {
		for (auto j : makeIntervalRange(reg[Species::Complex])) {
			radius +=
				pow(0.75 * prefactor * latticeParameter * (double)j, 1.0 / 3.0);
		}
		return radius / reg[Species::Complex].length();
	}
	if (lo.isOnAxis(Species::Trap)) {
		return ::xolotl::core::fecrCoreRadius;
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
