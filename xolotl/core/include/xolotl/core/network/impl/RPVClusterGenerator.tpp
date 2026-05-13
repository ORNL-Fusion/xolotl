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

	// Count the number of axis it is on
	int nAxis = 0;
	for (auto s : NetworkType::getSpeciesRange()) {
		if (region[s].begin() > 0) {
			nAxis++;
		}
	}

	// Cannot be on more than 1 axis
	if (nAxis > 1) {
		for (auto& r : result) {
			r = false;
		}
		return false;
	}

	// Cannot be 0 axis
	if (nAxis == 0)
		return true;

	// Check if others begin at 0
	auto othersBeginAtZero = [](const Region& reg, Species species) {
		for (auto s : NetworkType::getSpeciesRange()) {
			if (s.value != species && reg[s].begin() != 0) {
				return false;
			}
		}
		return true;
	};

	using detail::toIndex;

	// Get the bounds
	Composition lo = region.getOrigin();
	Composition hi = region.getUpperLimitPoint();

	// Smaller than the minimum size for grouping
	if (lo[Species::V] < _groupingMin && lo[Species::I] < _groupingMin &&
		lo[Species::Loop] < _groupingMin) {
		return true;
	}

	// I is grouped on its own
	if (lo[Species::I] > 0) {
		if (lo[Species::I] < _groupingMin &&
			othersBeginAtZero(region, Species::I)) {
			return true;
		}
		if (region[Species::I].end() > _maxSize) {
			return true;
		}
		if (region[Species::I].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::I].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::I)] = false;
		}
	}

	// V is grouped on its own
	if (lo[Species::V] > 0) {
		if (lo[Species::V] < _groupingMin &&
			othersBeginAtZero(region, Species::V)) {
			return true;
		}
		if (region[Species::V].end() > _maxSize) {
			return true;
		}
		if (region[Species::V].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::V].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::V)] = false;
		}
	}

	// Loop is grouped on its own
	if (lo[Species::Loop] > 0) {
		if (lo[Species::Loop] < _groupingMin &&
			othersBeginAtZero(region, Species::Loop)) {
			return true;
		}
		if (region[Species::Loop].end() > _maxSize) {
			return true;
		}
		if (region[Species::Loop].length() <
			util::max((double)(_groupingWidth + 1),
				pow(region[Species::Loop].begin(), 0.75) * 0.5)) {
			result[toIndex(Species::Loop)] = false;
		}
	}

	// Edges
	if (region[Species::I].end() > _maxSize + 1) {
		return true;
	}
	if (region[Species::V].end() > _maxSize + 1) {
		return true;
	}
	if (region[Species::Loop].end() > _maxSize + 1) {
		return true;
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
RPVClusterGenerator::select(const Region& region) const
{
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
	if (region[Species::I].begin() > _maxSize) {
		return false;
	}

	// Vacancies
	if (region[Species::V].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::V)) {
		return false;
	}
	if (region[Species::V].begin() > _maxSize) {
		return false;
	}

	// Loops
	if (region[Species::Loop].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::Loop)) {
		return false;
	}
	if (region[Species::Loop].begin() > 0 && region.isSimplex() &&
		region[Species::Loop].begin() <= _maxI) {
		return false;
	}
	if (region[Species::Loop].begin() > _maxSize) {
		return false;
	}

	if (region[Species::V].begin() == 0 && region[Species::I].begin() == 0 &&
		region[Species::Loop].end() - 1 <= _maxI &&
		region[Species::Loop].begin() > 1)
		return false;

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
RPVClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// V migration energies in eV
	constexpr Kokkos::Array<double, 5> vMigration = {
		0.0, 0.63, 0.50, 0.36, 0.44};

	const auto& reg = cluster.getRegion();
	Composition comp(reg.getOrigin());
	double migrationEnergy = util::infinity<double>;
	if (reg.isSimplex()) {
		if (comp.isOnAxis(Species::I)) {
			switch (comp[Species::I]) {
			case 1:
				return 0.34;
			case 2:
				return 0.42;
			case 3:
				return 0.43;
			default:
				return 0.9;
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto amtV = comp[Species::V];
			if (amtV < vMigration.size()) {
				migrationEnergy = vMigration[amtV];
			}
		}
	}
	else {
		if (comp.isOnAxis(Species::I))
			migrationEnergy = 0.9;
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
	Composition comp(reg.getOrigin());
	double diffusionFactor = 0.0;
	if (reg.isSimplex()) {
		if (comp.isOnAxis(Species::I)) {
			if (comp[Species::I] < 4) {
				diffusionFactor = iOneDiffusionFactor / comp[Species::I];
			}
			else {
				diffusionFactor = 2.965e11 * pow(comp[Species::I], -0.7);
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto amtV = comp[Species::V];
			if (amtV < vDiffusion.size()) {
				diffusionFactor = vDiffusion[amtV];
			}
		}
	}
	else {
		if (comp.isOnAxis(Species::I)) {
			Composition hiComp(reg.getUpperLimitPoint());
			double nF =
				(double)(comp[Species::I] + hiComp[Species::I] - 1) / 2.0;
			diffusionFactor = 2.965e11 * pow(comp[Species::I], -0.7);
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
	Composition comp(reg.getOrigin());
	Composition hiComp(reg.getUpperLimitPoint());

	// Constants
	const double prefactor =
		0.5 * latticeParameter * latticeParameter / ::xolotl::core::pi;
	constexpr double fecrBurgers = 0.8660254038;
	constexpr double fecrLoopBurgers = 1.0;

	// I case
	if (comp.isOnAxis(Species::I)) {
		// Sphere
		if (comp[Species::I] < 4) {
			radius = latticeParameter *
				cbrt(3.0 * comp[Species::I] / ::xolotl::core::pi) * 0.5;
		}
		// (111) Loop
		else {
			for (auto j : makeIntervalRange(reg[Species::I])) {
				radius += sqrt(((double)j * prefactor) / fecrBurgers);
			}
			// Average the radius
			radius /= reg[Species::I].length();
		}
	}
	// (100) Loop
	else if (comp.isOnAxis(Species::Loop)) {
		for (auto j : makeIntervalRange(reg[Species::Loop])) {
			radius += sqrt(((double)j * prefactor) / fecrLoopBurgers);
		}
		radius /= reg[Species::Loop].length();
	}
	// V case
	else {
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
