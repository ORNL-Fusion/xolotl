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
FeClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	result[0] = true;
	result[1] = true;
	result[2] = true;

	// I is never grouped
	if (region[Species::I].begin() > 0) {
		return true;
	}

	// V is never grouped
	if (region[Species::V].end() > 1 && region[Species::V].begin() < 11 &&
		region[Species::He].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// He is never grouped
	if (region[Species::He].end() > 1 && region[Species::He].begin() < 9 &&
		region[Species::V].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// HeV
	if (region[Species::He].begin() < _groupingMin &&
		region[Species::V].begin() < _groupingMin) {
		return true;
	}

	// Edges
	if (region[Species::He].end() > _maxHe + 1 &&
		region[Species::V].end() > _maxV + 1) {
		return true;
	}

	// Middle
	Composition lo = region.getOrigin();
	Composition hi = region.getUpperLimitPoint();
	auto amtHe = 0.5 * (lo[Species::He] + hi[Species::He] - 1);
	auto amtV = 0.5 * (lo[Species::V] + hi[Species::V] - 1);
	double amt = sqrt(amtHe * amtHe + amtV * amtV);
	double ibe = 4.88 +
		2.59 * (cbrt(amtV * amtV) - cbrt((amtV - 1.0) * (amtV - 1.0))) -
		2.5 * log(1.0 + (amtHe / amtV));
	auto distance = fabs(ibe - 1.5);
	//	if (distance < 1.2) {
	if (distance < 1.5) {
		auto comp = amt * amt * amt * 1.0e-6;
		if (region[Species::He].length() <
			util::max(_groupingWidthHe + 1.0, comp)) {
			result[0] = false;
		}
		if (region[Species::V].length() <
			util::max(_groupingWidthV + 1.0, comp)) {
			result[1] = false;
		}
	}
	else {
		auto comp = amt * amt * amt * 1.0e-4;
		if (region[Species::He].length() <
			util::max(_groupingWidthHe + 1.0, comp)) {
			result[0] = false;
		}
		if (region[Species::V].length() <
			util::max(_groupingWidthV + 1.0, comp)) {
			result[1] = false;
		}
	}

	// Edges
	if (region[Species::He].begin() == 0) {
		result[0] = true;
	}
	if (region[Species::V].begin() == 0) {
		result[1] = true;
	}
	if (region[Species::He].end() > _maxHe + 1) {
		result[0] = true;
	}
	if (region[Species::V].end() > _maxV + 1) {
		result[1] = true;
	}

	if (!result[0] && !result[1]) {
		return false;
	}

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
FeClusterGenerator::select(const Region& region) const
{
	// Remove 0
	if (region[Species::He].end() == 1 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}

	// Interstitials
	if (region[Species::I].begin() > 0 &&
		(region[Species::He].begin() > 0 || region[Species::V].begin() > 0)) {
		return false;
	}

	// Helium
	if (region[Species::He].begin() > 8 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}
	if (region[Species::He].begin() > _maxHe) {
		return false;
	}

	// Vacancy
	if (region[Species::V].begin() > 10 && region[Species::He].end() == 1 &&
		region[Species::I].end() == 1) {
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
FeClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// I migration energy in eV
	constexpr double iOneMigrationEnergy = 0.34;
	// He migration energies in eV
	constexpr Kokkos::Array<double, 4> heMigration = {0.0, 0.06, 0.06, 0.06};
	// V migration energies in eV
	constexpr Kokkos::Array<double, 5> vMigration = {
		0.0, 0.67, 0.62, 0.37, 0.48};

	const auto& reg = cluster.getRegion();
	double migrationEnergy = util::infinity<double>;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			if (comp[Species::I] == 1) {
				migrationEnergy = iOneMigrationEnergy;
			}
		}
		else if (comp.isOnAxis(Species::He)) {
			auto amtHe = comp[Species::He];
			if (amtHe < heMigration.size()) {
				migrationEnergy = heMigration[amtHe];
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
FeClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	// I diffusion factors in nm^2/s
	constexpr double iOneDiffusionFactor = 1.0e+11;
	// He diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 4> heDiffusion = {
		0.0, 1.0e+11, 5.0e+10, 3.3e+10};
	// V diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 5> vDiffusion = {
		0.0, 1.0e+11, 5.0e+10, 3.3e+10, 2.5e+10};

	const auto& reg = cluster.getRegion();
	double diffusionFactor = 0.0;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			if (comp[Species::I] == 1) {
				diffusionFactor = iOneDiffusionFactor;
			}
		}
		else if (comp.isOnAxis(Species::He)) {
			auto amtHe = comp[Species::He];
			if (amtHe < heDiffusion.size()) {
				diffusionFactor = heDiffusion[amtHe];
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
FeClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
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
		else if (comp.isOnAxis(Species::He)) {
			double FourPi = 4.0 * ::xolotl::core::pi;
			double aCubed = pow(latticeParameter, 3);
			double termOne =
				pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::He],
					(1.0 / 3.0));
			double termTwo =
				pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
			radius = impurityRadius + termOne - termTwo;
		}
		else {
			radius = latticeParameter *
				pow((3.0 * comp[Species::V]) / ::xolotl::core::pi,
					(1.0 / 3.0)) *
				0.5;
		}
	}
	else {
		// Loop on the V range
		for (auto j : makeIntervalRange(reg[Species::V])) {
			radius += latticeParameter *
				pow((3.0 * (double)j) / ::xolotl::core::pi, (1.0 / 3.0)) * 0.5;
		}
		// Average the radius
		radius /= reg[Species::V].length();
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
