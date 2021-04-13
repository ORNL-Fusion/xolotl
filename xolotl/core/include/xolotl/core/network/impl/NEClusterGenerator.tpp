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
	if (region[Species::V].end() > 1 && region[Species::V].begin() < 7 &&
		region[Species::Xe].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// Xe is never grouped
	if (region[Species::Xe].end() > 1 && region[Species::Xe].begin() < 2 &&
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
	if (region[Species::Xe].length() <
		util::max((double)(_groupingWidthXe + 1),
			region[Species::Xe].begin() * 1.0e-1)) {
		result[0] = false;
	}
	if (region[Species::V].length() <
		util::max((double)(_groupingWidthV + 1),
			region[Species::V].begin() * 1.0e-1)) {
		result[1] = false;
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
	if (region[Species::Xe].begin() > 1 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}
	if (region[Species::Xe].begin() > _maxXe) {
		return false;
	}

	// Vacancy
	if (region[Species::V].begin() > 6 && region[Species::Xe].end() == 1 &&
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
	// Not used?
	return 0.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
NEClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// I migration energy in eV
	constexpr double iOneMigrationEnergy = 3.5;
	// Xe migration energies in eV
	constexpr double xeOneMigrationEnergy = 1.0;
	// V migration energies in eV
	constexpr double vOneMigrationEnergy = 3.5;
	// XeV migration energies in eV
	constexpr Kokkos::Array<double, 9> xevMigrationEnergies = {
		util::infinity<double>, util::infinity<double>, 3.70, 4.99, 2.77, 7.37,
		10.13, 5.43, 3.13};

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
			if (comp[Species::Xe] == 1) {
				migrationEnergy = xeOneMigrationEnergy;
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			if (comp[Species::V] == 1) {
				migrationEnergy = vOneMigrationEnergy;
			}
		}
		else if (comp[Species::V] < xevMigrationEnergies.size() &&
			comp[Species::Xe] == 1) {
			migrationEnergy = xevMigrationEnergies[comp[Species::V]];
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
	constexpr double iOneDiffusionFactor = 4.077e+11;
	// Xe migration energies in eV
	constexpr double xeOneDiffusionFactor = 1.247e+10;
	// V migration energies in eV
	constexpr double vOneDiffusionFactor = 4.249e+11;
	// XeV migration energies in eV
	constexpr Kokkos::Array<double, 9> xevDiffusionFactors = {0.0, 0.0,
		6.483e+11, 7.48e+10, 9.974e+10, 2.493e+10, 1.995e+11, 2.493e+10,
		2.493e+10};

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
			if (comp[Species::Xe] == 1) {
				diffusionFactor = xeOneDiffusionFactor;
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			if (comp[Species::V] == 1) {
				diffusionFactor = vOneDiffusionFactor;
			}
		}
		else if (comp[Species::V] < xevDiffusionFactors.size() &&
			comp[Species::Xe] == 1) {
			diffusionFactor = xevDiffusionFactors[comp[Species::V]];
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
	double FourPi = 4.0 * ::xolotl::core::pi;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			radius = latticeParameter * cbrt(1.0 / ::xolotl::core::pi) * 0.5;
		}
		else if (comp.isOnAxis(Species::Xe)) {
			radius = impurityRadius;
		}
		else if (comp.isOnAxis(Species::V)) {
			radius = latticeParameter *
				cbrt((comp[Species::V]) / ::xolotl::core::pi) * 0.5;
		}
		else {
			radius =
				pow((3.0 * (double)comp[Species::Xe]) / (FourPi * _density),
					(1.0 / 3.0));
		}
	}
	else {
		// Loop on the Xe range
		for (auto j : makeIntervalRange(reg[Species::Xe])) {
			radius += pow((3.0 * (double)j) / (FourPi * _density), (1.0 / 3.0));
		}
		// Average the radius
		radius /= reg[Species::Xe].length();
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
