#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
NEClusterGenerator::NEClusterGenerator(const options::IOptions& options) :
	_maxXe(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthXe(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB()),
	_xeDiffusivity(options.getXenonDiffusivity()),
	_xeDiffusive(_xeDiffusivity > 0.0),
	_density(options.getDensity())
{
}

NEClusterGenerator::NEClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxXe(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthXe(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB()),
	_xeDiffusivity(options.getXenonDiffusivity()),
	_xeDiffusive(_xeDiffusivity > 0.0),
	_density(options.getDensity())
{
}

KOKKOS_INLINE_FUNCTION
bool
NEClusterGenerator::refine(const Region& region, BoolArray& result) const
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
		region[Species::Xe].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// Xe is never grouped
	if (region[Species::Xe].end() > 1 && region[Species::Xe].begin() < 9 &&
		region[Species::V].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// XeV
	if (region[Species::Xe].begin() < _groupingMin &&
		region[Species::V].begin() < _groupingMin) {
		return true;
	}

	// Edges
	if (region[Species::Xe].end() > _maxXe + 1 &&
		region[Species::V].end() > _maxV + 1) {
		return true;
	}

	// Middle
	Composition lo = region.getOrigin();
	Composition hi = region.getUpperLimitPoint();
	auto amtXe = 0.5 * (lo[Species::Xe] + hi[Species::Xe] - 1);
	auto amtV = 0.5 * (lo[Species::V] + hi[Species::V] - 1);
	double amt = sqrt(amtXe * amtXe + amtV * amtV);
	double ibe = 4.88 +
		2.59 * (cbrt(amtV * amtV) - cbrt((amtV - 1.0) * (amtV - 1.0))) -
		2.5 * log(1.0 + (amtXe / amtV));
	auto distance = fabs(ibe - 1.5);
	//	if (distance < 1.2) {
	if (distance < 1.5) {
		auto comp = amt * amt * amt * 1.0e-6;
		if (region[Species::Xe].length() <
			util::max(_groupingWidthXe + 1.0, comp)) {
			result[0] = false;
		}
		if (region[Species::V].length() <
			util::max(_groupingWidthV + 1.0, comp)) {
			result[1] = false;
		}
	}
	else {
		auto comp = amt * amt * amt * 1.0e-4;
		if (region[Species::Xe].length() <
			util::max(_groupingWidthXe + 1.0, comp)) {
			result[0] = false;
		}
		if (region[Species::V].length() <
			util::max(_groupingWidthV + 1.0, comp)) {
			result[1] = false;
		}
	}

	// Edges
	if (region[Species::Xe].begin() == 0) {
		result[0] = true;
	}
	if (region[Species::V].begin() == 0) {
		result[1] = true;
	}
	if (region[Species::Xe].end() > _maxXe + 1) {
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
NEClusterGenerator::select(const Region& region) const
{
	// Remove 0
	if (region[Species::Xe].end() == 1 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}

	// Interstitials
	if (region[Species::I].begin() > 0 &&
		(region[Species::Xe].begin() > 0 || region[Species::V].begin() > 0)) {
		return false;
	}

	// Xenon
	if (region[Species::Xe].begin() > 8 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}
	if (region[Species::Xe].begin() > _maxXe) {
		return false;
	}

	// Vacancy
	if (region[Species::V].begin() > 10 && region[Species::Xe].end() == 1 &&
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
NEClusterGenerator::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	/**
	 * The set of xenon formation energies up to Xe_29 indexed by size. That is
	 * E_(f,Xe_1) = xeFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	constexpr Kokkos::Array<double, 30> xeFormationEnergies = {0.0, 7.0, 12.15,
		17.15, 21.90, 26.50, 31.05, 35.30, 39.45, 43.00, 46.90, 50.65, 53.90,
		56.90, 59.80, 62.55, 65.05, 67.45, 69.45, 71.20, 72.75, 74.15, 75.35,
		76.40, 77.25, 77.95, 78.45, 78.80, 78.95, 79.0};

	const auto& reg = cluster.getRegion();
	if (reg.isSimplex()) {
		auto amtXe = reg.getOrigin()[0];
		if (amtXe < xeFormationEnergies.size()) {
			return xeFormationEnergies[amtXe];
		}
		else {
			return 79.0;
		}
	}
	return 79.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
NEClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// I migration energy in eV
	constexpr double iOneMigrationEnergy = 0.34;
	// Xe migration energies in eV
	constexpr Kokkos::Array<double, 4> xeMigration = {0.0, 0.06, 0.06, 0.06};
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
		else if (comp.isOnAxis(Species::Xe)) {
			auto amtXe = comp[Species::Xe];
			if (amtXe < xeMigration.size()) {
				migrationEnergy = xeMigration[amtXe];
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
NEClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	// I diffusion factors in nm^2/s
	constexpr double iOneDiffusionFactor = 1.0e+11;
	// Xe diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 4> xeDiffusion = {
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
		else if (comp.isOnAxis(Species::Xe)) {
			auto amtXe = comp[Species::Xe];
			if (amtXe < xeDiffusion.size()) {
				diffusionFactor = xeDiffusion[amtXe];
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
NEClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
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
		else if (comp.isOnAxis(Species::Xe)) {
			double FourPi = 4.0 * ::xolotl::core::pi;
			double aCubed = pow(latticeParameter, 3);
			double termOne =
				pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::Xe],
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
