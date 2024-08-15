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
T91ClusterGenerator::refine(const Region& region, BoolArray& result) const
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
	if (region[Species::He].end() > 1 && region[Species::He].begin() < 5 &&
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
	double amtHe = 0.5 * (lo[Species::He] + hi[Species::He] - 1);
	double amtV = 0.5 * (lo[Species::V] + hi[Species::V] - 1);
	double amt = sqrt(amtHe * amtHe + amtV * amtV);
	double highRatio =
		(hi[Species::He] - 1) / xolotl::util::max(1.0, (double)lo[Species::V]);
	double lowRatio =
		lo[Species::He] / xolotl::util::max(1.0, (double)(hi[Species::V] - 1));
	if (lowRatio > 3.2 or highRatio < 1.8) {
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
	else {
		auto comp = amt * 1.0e-2;
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
T91ClusterGenerator::select(const Region& region) const
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
	if (region[Species::He].begin() > 4 && region[Species::V].end() == 1 &&
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
T91ClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// I migration energy in eV
	constexpr double iOneMigrationEnergy = 0.2;
	// He migration energies in eV
	constexpr Kokkos::Array<double, 5> heMigration = {
		0.0, 0.064, 0.064, 0.064, 0.054};
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
		else if (comp[Species::V] == 1) {
			if (comp[Species::He] == 2)
				migrationEnergy = 0.31;
			else if (comp[Species::He] == 3)
				migrationEnergy = 0.31;
			else if (comp[Species::He] == 4)
				migrationEnergy = 0.31;
		}
		else if (comp[Species::V] == 2 and comp[Species::He] == 1)
			migrationEnergy = 1.17;
		else if (comp[Species::V] == 3 and comp[Species::He] == 2)
			migrationEnergy = 0.53;
	}
	return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
T91ClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	// I diffusion factors in nm^2/s
	constexpr double iOneDiffusionFactor = 2.8e+10;
	// He diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 5> heDiffusion = {
		0.0, 5.2e+11, 3.1e+10, 8.0e+9, 1.2e+9};
	// V diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 5> vDiffusion = {
		0.0, 6.8e+13, 6.8e+13, 6.8e+13, 6.8e+13};

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
		else if (comp[Species::V] == 1) {
			if (comp[Species::He] == 2)
				diffusionFactor = 3.34e+11;
			else if (comp[Species::He] == 3)
				diffusionFactor = 3.2e+10;
			else if (comp[Species::He] == 4)
				diffusionFactor = 2.13e+9;
		}
		else if (comp[Species::V] == 2 and comp[Species::He] == 1)
			diffusionFactor = 9.93e+9;
		else if (comp[Species::V] == 3 and comp[Species::He] == 2)
			diffusionFactor = 9.93e+9;
	}

	return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
T91ClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
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
