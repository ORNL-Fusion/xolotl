#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TSpeciesEnum>
PSIClusterGenerator<TSpeciesEnum>::PSIClusterGenerator(
	const options::IOptions& opts) :
	_hydrogenRadiusFactor(opts.getHydrogenFactor()),
	_maxHe(opts.getMaxImpurity()),
	_maxD(opts.getMaxD()),
	_maxT(opts.getMaxT()),
	_maxV(opts.getMaxV()),
	_maxI(opts.getMaxI()),
	_groupingMin(opts.getGroupingMin()),
	_groupingWidthA(opts.getGroupingWidthA()),
	_groupingWidthB(opts.getGroupingWidthB()),
	_hevRatio(opts.getHeVRatio())
{
}

template <typename TSpeciesEnum>
PSIClusterGenerator<TSpeciesEnum>::PSIClusterGenerator(
	const options::IOptions& opts, std::size_t refineDepth) :
	Superclass(refineDepth),
	_hydrogenRadiusFactor(opts.getHydrogenFactor()),
	_maxHe(opts.getMaxImpurity()),
	_maxD(opts.getMaxD()),
	_maxT(opts.getMaxT()),
	_maxV(opts.getMaxV()),
	_maxI(opts.getMaxI()),
	_groupingMin(opts.getGroupingMin()),
	_groupingWidthA(opts.getGroupingWidthA()),
	_groupingWidthB(opts.getGroupingWidthB()),
	_hevRatio(opts.getHeVRatio())
{
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
bool
PSIClusterGenerator<TSpeciesEnum>::refine(
	const Region& region, BoolArray& result) const
{
	using detail::toIndex;
	using psi::getMaxHePerV;
	using psi::hasDeuterium;
	using psi::hasTritium;

	for (auto& r : result) {
		r = true;
	}

	Composition lo = region.getOrigin();
	Composition hi = region.getUpperLimitPoint();

	auto othersBeginAtZero = [](const Region& reg, Species species) {
		for (auto s : NetworkType::getSpeciesRange()) {
			if (s.value != species && reg[s].begin() != 0) {
				return false;
			}
		}
		return true;
	};

	if (lo[Species::I] > 0) {
		if (lo[Species::I] < _groupingMin &&
			othersBeginAtZero(region, Species::I)) {
			return true;
		}
		if (region[Species::I].end() > _maxI) {
			return true;
		}
		if (region[Species::I].length() <
			util::max((double)(_groupingWidthB + 1),
				region[Species::I].begin() * 1.0e-2)) {
			result[toIndex(Species::I)] = false;
			return true;
		}
	}

	// No impurity case
	if (lo[Species::V] > 0 and _maxHe == 0 and _maxD == 0 and _maxT == 0) {
		if (lo[Species::V] < _groupingMin &&
			othersBeginAtZero(region, Species::V)) {
			return true;
		}
		if (region[Species::V].end() > _maxV) {
			return true;
		}
		if (region[Species::V].length() <
			util::max((double)(_groupingWidthB + 1),
				region[Species::V].begin() * 1.0e-2)) {
			result[toIndex(Species::V)] = false;
			return true;
		}
	}
	else if (_maxV < 10000) {
		// V is never grouped
		if (hi[Species::V] > 1 && othersBeginAtZero(region, Species::V)) {
			return true;
		}
	}

	// Don't group under the given min for V
	if (lo[Species::V] < _groupingMin) {
		return true;
	}

	// Skip the rest if there is no impurities
	if (_maxHe == 0 and _maxD == 0 and _maxT == 0)
		return true;

	// He is never grouped
	if (hi[Species::He] > 1 && othersBeginAtZero(region, Species::He)) {
		return true;
	}

	// D is never grouped
	if constexpr (hasDeuterium<Species>) {
		if (hi[Species::D] > 1 && othersBeginAtZero(region, Species::D)) {
			return true;
		}
	}

	// T is never grouped
	if constexpr (hasTritium<Species>) {
		if (hi[Species::T] > 1 && othersBeginAtZero(region, Species::T)) {
			return true;
		}
	}

	// Else refine around the edge
	auto maxDPerV = [hevRatio = _hevRatio, maxD = _maxD](AmountType amtV) {
		return (2.0 / 3.0) * getMaxHePerV(amtV, hevRatio) * (maxD > 0);
	};
	auto maxTPerV = [hevRatio = _hevRatio, maxT = _maxT](AmountType amtV) {
		return (2.0 / 3.0) * getMaxHePerV(amtV, hevRatio) * (maxT > 0);
	};
	if (hi[Species::V] > 1) {
		double factor = 1.0e-1;

		if (lo[Species::He] <= getMaxHePerV(hi[Species::V] - 1, _hevRatio) &&
			hi[Species::He] - 1 >=
				getMaxHePerV(lo[Species::V] - 1, _hevRatio)) {
			if constexpr (hasDeuterium<Species>) {
				if (region[Species::D].length() <
					util::max((double)(_groupingWidthA + 1),
						(double)lo[Species::D] * (double)lo[Species::D] *
							factor)) {
					result[toIndex(Species::D)] = false;
				}
				if (lo[Species::D] <= maxDPerV(hi[Species::V] - 1) &&
					hi[Species::D] - 1 >= maxDPerV(lo[Species::V] - 1) &&
					_maxD > 0) {
					result[toIndex(Species::D)] = true;
				}
			}

			if constexpr (hasTritium<Species>) {
				if (region[Species::T].length() <
					util::max((double)(_groupingWidthA + 1),
						(double)lo[Species::T] * (double)lo[Species::T] *
							factor)) {
					result[toIndex(Species::T)] = false;
				}
				if (lo[Species::T] <= maxTPerV(hi[Species::V] - 1) &&
					hi[Species::T] - 1 >= maxTPerV(lo[Species::V] - 1) &&
					_maxT > 0) {
					result[toIndex(Species::T)] = true;
				}
			}

			if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
				auto hiH = hi[Species::D] + hi[Species::T] - 2;
				if (hiH > (2.0 / 3.0) * lo[Species::He] + 0.5) {
					result[toIndex(Species::D)] = true;
					result[toIndex(Species::T)] = true;
				}
			}

			// Border case
			if constexpr (hasDeuterium<Species>) {
				if (lo[Species::D] == 0) {
					result[toIndex(Species::D)] = true;
				}
			}
			if constexpr (hasTritium<Species>) {
				if (lo[Species::T] == 0) {
					result[toIndex(Species::T)] = true;
				}
			}

			// Doesn't fully refine only for large networks
			if (_maxV < 10000)
				return true;

			if (region[Species::He].length() <=
				(lo[Species::He] - 4.0 * _groupingMin) * factor)
				result[toIndex(Species::He)] = false;
			if (hi[Species::He] > getMaxHePerV(_maxV, _hevRatio) + 1) {
				result[toIndex(Species::He)] = true;
			}

			if (region[Species::V].length() <=
				(lo[Species::V] - _groupingMin) * factor)
				result[toIndex(Species::V)] = false;

			if (hi[Species::V] > _maxV + 1) {
				result[toIndex(Species::V)] = true;
			}

			return true;
		}
	}

	double factor = 5.0e-1;

	if (_maxV < 10000) {
		if (region[Species::V].length() < _groupingWidthB + 1) {
			result[toIndex(Species::V)] = false;
		}
	}
	else {
		if (region[Species::V].length() <
			util::max((double)(_groupingWidthB + 1),
				(double)lo[Species::V] * (double)lo[Species::V] * factor)) {
			result[toIndex(Species::V)] = false;
		}
	}

	if (region[Species::He].length() <
		util::max((double)(_groupingWidthA + 1),
			(double)lo[Species::He] * (double)lo[Species::He] * factor)) {
		result[toIndex(Species::He)] = false;
	}

	if constexpr (hasDeuterium<Species>) {
		if (region[Species::D].length() <
			util::max((double)(_groupingWidthA + 1),
				(double)lo[Species::D] * (double)lo[Species::D] * factor)) {
			result[toIndex(Species::D)] = false;
		}
	}

	if constexpr (hasTritium<Species>) {
		if (region[Species::T].length() <
			util::max((double)(_groupingWidthA + 1),
				(double)lo[Species::T] * (double)lo[Species::T] * factor)) {
			result[toIndex(Species::T)] = false;
		}
	}

	// No He
	if (lo[Species::He] == 0) {
		if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
			auto hiH = hi[Species::D] + hi[Species::T] - 2;
			if (hiH > 6 * (hi[Species::V] - 1)) {
				result[toIndex(Species::D)] = true;
				result[toIndex(Species::T)] = true;
			}
		}
		else if constexpr (hasDeuterium<Species>) {
			if (hi[Species::D] - 1 > 6 * (hi[Species::V] - 1)) {
				result[toIndex(Species::D)] = true;
			}
		}
		else if constexpr (hasTritium<Species>) {
			if (hi[Species::T] - 1 > 6 * (hi[Species::V] - 1)) {
				result[toIndex(Species::T)] = true;
			}
		}
	}

	// Border case
	if (hi[Species::V] > _maxV + 1) {
		result[toIndex(Species::V)] = true;
	}
	if (lo[Species::V] == 0) {
		result[toIndex(Species::V)] = true;
	}
	if (lo[Species::He] == 0) {
		result[toIndex(Species::He)] = true;
	}
	if constexpr (hasDeuterium<Species>) {
		if (lo[Species::D] == 0) {
			result[toIndex(Species::D)] = true;
		}
	}
	if constexpr (hasTritium<Species>) {
		if (lo[Species::T] == 0) {
			result[toIndex(Species::T)] = true;
		}
	}

	int axis = 0;
	for (auto& r : result) {
		axis += r;
	}

	if (axis == 0)
		return false;

	return true;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
bool
PSIClusterGenerator<TSpeciesEnum>::select(const Region& region) const
{
	using psi::getMaxHePerV;
	using psi::hasDeuterium;
	using psi::hasTritium;

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

	// Interstitials
	if (region[Species::I].begin() > 0 &&
		!region.getOrigin().isOnAxis(Species::I)) {
		return false;
	}
	if (region[Species::I].begin() > _maxI) {
		return false;
	}

	// Helium
	if (region[Species::He].begin() > _maxHe &&
		othersEndAtOne(region, Species::He)) {
		return false;
	}

	// Deuterium
	if constexpr (hasDeuterium<Species>) {
		if (region[Species::D].begin() > _maxD &&
			othersEndAtOne(region, Species::D)) {
			return false;
		}
	}

	// Tritium
	if constexpr (hasTritium<Species>) {
		if (region[Species::T].begin() > _maxT &&
			othersEndAtOne(region, Species::T)) {
			return false;
		}
	}

	// Vacancy
	if (region[Species::V].begin() > 0 and
		!region.getOrigin().isOnAxis(Species::V) and _maxHe == 0 and
		_maxD == 0 and _maxT == 0) {
		return false;
	}
	if (region[Species::V].begin() > _maxV) {
		return false;
	}

	// Can't cluster without V
	if (region[Species::V].end() == 1) {
		if (region[Species::He].begin() > _maxHe)
			return false;
		if constexpr (hasDeuterium<Species>) {
			if (region[Species::D].begin() > _maxD)
				return false;
			if (region[Species::He].begin() > 0 &&
				region[Species::D].begin() > 0) {
				return false;
			}
		}
		if constexpr (hasTritium<Species>) {
			if (region[Species::T].begin() > _maxT)
				return false;
			if (region[Species::He].begin() > 0 &&
				region[Species::T].begin() > 0) {
				return false;
			}
		}
		if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
			if (region[Species::D].begin() > 0 &&
				region[Species::T].begin() > 0) {
				return false;
			}
		}
	}

	auto maxHPerV = [hevRatio = _hevRatio](AmountType amtV) {
		return (2.0 / 3.0) * getMaxHePerV(amtV, hevRatio);
	};

	// The edge
	if (region[Species::V].end() > 1) {
		Composition lo = region.getOrigin();
		Composition hi = region.getUpperLimitPoint();
		auto hiV = util::min(hi[Species::V] - 1, _maxV);
		auto hiHe =
			util::min(hi[Species::He] - 1, getMaxHePerV(_maxV, _hevRatio));

		// Too many helium
		if (lo[Species::He] > getMaxHePerV(hiV, _hevRatio)) {
			return false;
		}

		// Too many deuterium
		if constexpr (hasDeuterium<Species>) {
			if (lo[Species::D] > maxHPerV(hiV)) {
				return false;
			}
		}

		// Too many tritium
		if constexpr (hasTritium<Species>) {
			if (lo[Species::T] > maxHPerV(hiV)) {
				return false;
			}
		}

		// Too many hydrogen
		if constexpr (hasDeuterium<Species> && hasTritium<Species>) {
			auto loH = lo[Species::D] + lo[Species::T];
			if (lo[Species::He] == 0) {
				if (loH > 6 * hiV) {
					return false;
				}
			}
			else {
				if (loH > (2.0 / 3.0) * hiHe + 0.5) {
					return false;
				}
			}
		}
		else if constexpr (hasDeuterium<Species>) {
			if (lo[Species::He] == 0) {
				if (lo[Species::D] > 6 * hiV) {
					return false;
				}
			}
			else {
				if (lo[Species::D] > (2.0 / 3.0) * hiHe + 0.5) {
					return false;
				}
			}
		}
		else if constexpr (hasTritium<Species>) {
			if (lo[Species::He] == 0) {
				if (lo[Species::T] > 6 * hiV) {
					return false;
				}
			}
			else {
				if (lo[Species::T] > (2.0 / 3.0) * hiHe + 0.5) {
					return false;
				}
			}
		}
	}

	return true;
}

template <typename TSpeciesEnum>
template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<TSpeciesEnum>::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	constexpr auto infinity = util::infinity<double>;

	// I formation energies in eV
	constexpr Kokkos::Array<double, 7> iFormationEnergies = {
		0.0, 10.0, 18.5, 27.0, 35.0, 42.5, 48.0};
	// He formation energies in eV
	constexpr Kokkos::Array<double, 9> heFormationEnergies = {
		0.0, 6.15, 11.44, 16.35, 21.0, 26.1, 30.24, 34.93, 38.80};

	const auto& reg = cluster.getRegion();
	double formationEnergy{};
	// All the possibly grouped ones are 0
	if (reg[Species::V].end() >= _groupingMin)
		return formationEnergy;

	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			auto amtI = comp[Species::I];
			if (amtI < iFormationEnergies.size()) {
				return iFormationEnergies[amtI];
			}
			else {
				return iFormationEnergies[iFormationEnergies.size() - 1] +
					(double)(amtI + 1 - iFormationEnergies.size());
			}
		}
		if (comp.isOnAxis(Species::He)) {
			auto amtHe = comp[Species::He];
			if (amtHe < heFormationEnergies.size()) {
				return heFormationEnergies[amtHe];
			}
			else {
				return infinity;
			}
		}
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp.isOnAxis(Species::D)) {
				return infinity;
			}
		}
		if constexpr (psi::hasTritium<Species>) {
			if (comp.isOnAxis(Species::T)) {
				return infinity;
			}
		}

		return getHeVFormationEnergy(comp);
	}

	return formationEnergy;
}

template <typename TSpeciesEnum>
template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<TSpeciesEnum>::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// I migration energies in eV
	constexpr Kokkos::Array<double, 6> iMigration = {
		0.0, 0.01, 0.02, 0.03, 0.04, 0.05};
	// He migration energies in eV
	constexpr Kokkos::Array<double, 8> heMigration = {
		0.0, 0.13, 0.20, 0.25, 0.20, 0.12, 0.3, 0.4};
	// The migration energy for a single deuterium.
	constexpr double dOneMigrationEnergy = 0.38;
	// The migration energy for a single tritium.
	static constexpr double tOneMigrationEnergy = 0.38;
	// The migration energy for a single vacancy in eV
	static constexpr double vOneMigration = 1.30;

	const auto& reg = cluster.getRegion();
	double migrationEnergy = util::infinity<double>;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			auto amtI = comp[Species::I];
			if (amtI < iMigration.size()) {
				return iMigration[amtI];
			}
			else {
				return util::min((double)amtI, 15.0) * 0.1;
			}
		}
		if (comp.isOnAxis(Species::He)) {
			auto amtHe = comp[Species::He];
			if (amtHe < heMigration.size()) {
				return heMigration[amtHe];
			}
		}
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp.isOnAxis(Species::D)) {
				if (comp[Species::D] == 1) {
					return dOneMigrationEnergy;
				}
			}
		}
		if constexpr (psi::hasTritium<Species>) {
			if (comp.isOnAxis(Species::T)) {
				if (comp[Species::T] == 1) {
					return tOneMigrationEnergy;
				}
			}
		}
		if (comp.isOnAxis(Species::V)) {
			if (comp[Species::V] <= 1) {
				return vOneMigration;
			}
		}
	}
	else {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			migrationEnergy = 0.0;
			for (auto j : makeIntervalRange(reg[Species::I])) {
				migrationEnergy += util::min((double)j, 15.0) * 0.1;
			}
			// Average the energy
			migrationEnergy /= reg[Species::I].length();
		}
	}

	return migrationEnergy;
}

template <typename TSpeciesEnum>
template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<TSpeciesEnum>::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	// I diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 6> iDiffusion = {
		0.0, 8.8e+10, 8.0e+10, 3.9e+10, 2.0e+10, 1.0e+10};
	// He diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 8> heDiffusion = {
		0.0, 2.9e+10, 3.2e+10, 2.3e+10, 1.7e+10, 5.0e+09, 1.0e+09, 5.0e+08};
	// The diffusion factor for a single deuterium.
	constexpr double dOneDiffusionFactor = 2.83e+11;
	// The diffusion factor for a single tritium.
	constexpr double tOneDiffusionFactor = 2.31e+11;
	// The diffusion factor for a single vacancy in nm^2/s
	constexpr double vOneDiffusion = 1.8e+12;

	const auto& reg = cluster.getRegion();
	double diffusionFactor = 0.0;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			auto amtI = comp[Species::I];
			if (amtI < iDiffusion.size()) {
				return iDiffusion[amtI];
			}
			else {
				return iDiffusion[1] / (double)amtI;
			}
		}
		if (comp.isOnAxis(Species::He)) {
			auto amtHe = comp[Species::He];
			if (amtHe < heDiffusion.size()) {
				return heDiffusion[amtHe];
			}
		}
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp.isOnAxis(Species::D)) {
				if (comp[Species::D] == 1) {
					return dOneDiffusionFactor;
				}
			}
		}
		if constexpr (psi::hasTritium<Species>) {
			if (comp.isOnAxis(Species::T)) {
				if (comp[Species::T] == 1) {
					return tOneDiffusionFactor;
				}
			}
		}
		if (comp.isOnAxis(Species::V)) {
			if (comp[Species::V] == 1) {
				return vOneDiffusion;
			}
		}
	}
	else {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			for (auto j : makeIntervalRange(reg[Species::I])) {
				diffusionFactor += iDiffusion[1] / (double)j;
			}
			// Average the energy
			diffusionFactor /= reg[Species::I].length();
		}
	}

	return diffusionFactor;
}

template <typename TSpeciesEnum>
template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<TSpeciesEnum>::getReactionRadius(
	const Cluster<PlsmContext>& cluster, double latticeParameter,
	double interstitialBias, double impurityRadius) const noexcept
{
	const auto& reg = cluster.getRegion();
	double radius = 0.0;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			double EightPi = 8.0 * ::xolotl::core::pi;
			double aCubed = pow(latticeParameter, 3.0);
			double termOne =
				interstitialBias * (sqrt(3.0) / 4.0) * latticeParameter;
			double termTwo =
				pow((3.0 / EightPi) * aCubed * comp[Species::I], (1.0 / 3.0));
			double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));
			return termOne + termTwo - termThree;
		}
		if (comp.isOnAxis(Species::He)) {
			double FourPi = 4.0 * ::xolotl::core::pi;
			double aCubed = pow(latticeParameter, 3);
			double termOne =
				pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::He],
					(1.0 / 3.0));
			double termTwo =
				pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
			return impurityRadius + termOne - termTwo;
		}
		if constexpr (psi::hasDeuterium<Species>) {
			if (comp.isOnAxis(Species::D)) {
				double FourPi = 4.0 * ::xolotl::core::pi;
				double aCubed = pow(latticeParameter, 3);
				double termOne = pow(
					(3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::D],
					(1.0 / 3.0));
				double termTwo =
					pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
				return (impurityRadius + termOne - termTwo) *
					_hydrogenRadiusFactor;
			}
		}
		if constexpr (psi::hasTritium<Species>) {
			if (comp.isOnAxis(Species::T)) {
				double FourPi = 4.0 * ::xolotl::core::pi;
				double aCubed = pow(latticeParameter, 3);
				double termOne = pow(
					(3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::T],
					(1.0 / 3.0));
				double termTwo =
					pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
				return (impurityRadius + termOne - termTwo) *
					_hydrogenRadiusFactor;
			}
		}

		return (sqrt(3.0) / 4.0) * latticeParameter +
			pow((3.0 * pow(latticeParameter, 3.0) * comp[Species::V]) /
					(8.0 * ::xolotl::core::pi),
				(1.0 / 3.0)) -
			pow((3.0 * pow(latticeParameter, 3.0)) / (8.0 * ::xolotl::core::pi),
				(1.0 / 3.0));
	}
	else {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			double EightPi = 8.0 * ::xolotl::core::pi;
			double aCubed = pow(latticeParameter, 3.0);
			double termOne =
				interstitialBias * (sqrt(3.0) / 4.0) * latticeParameter;
			double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));
			for (auto j : makeIntervalRange(reg[Species::I])) {
				double termTwo =
					pow((3.0 / EightPi) * aCubed * (double)j, (1.0 / 3.0));
				radius += termOne + termTwo - termThree;
			}
			// Average the energy
			radius /= reg[Species::I].length();
		}
		// Loop on the V range
		else {
			// Loop on the V range
			for (auto j : makeIntervalRange(reg[Species::V])) {
				radius += (sqrt(3.0) / 4.0) * latticeParameter +
					pow((3.0 * pow(latticeParameter, 3.0) * (double)j) /
							(8.0 * ::xolotl::core::pi),
						(1.0 / 3.0)) -
					pow((3.0 * pow(latticeParameter, 3.0)) /
							(8.0 * ::xolotl::core::pi),
						(1.0 / 3.0));
			}
			// Average the radius
			radius /= reg[Species::V].length();
		}
	}

	return radius;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<TSpeciesEnum>::getHeVFormationEnergy(Composition comp)
{
	// V formation energies in eV
	constexpr Kokkos::Array<double, 2> vFormationEnergies = {3.6, 7.25};

	// Coefficients for the Legendre polynomial fit
	// Low means V <= 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c0CoefficientsLow = {
		253.35, 435.36, 336.50, 198.92, 95.154, 21.544};
	// Coefficients for c_1 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c1CoefficientsLow = {
		493.29, 1061.3, 1023.9, 662.92, 294.24, 66.962};
	// Coefficients for c_2 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c2CoefficientsLow = {
		410.40, 994.89, 1044.6, 689.41, 286.52, 60.712};
	// Coefficients for c_3 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c3CoefficientsLow = {
		152.99, 353.16, 356.10, 225.75, 87.077, 15.640};
	// High means V > 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c0CoefficientsHigh = {
		-847.90, -3346.9, -4510.3, -3094.7, -971.18, -83.770};
	// Coefficients for c_1 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c1CoefficientsHigh = {
		-1589.3, -4894.6, -6001.8, -4057.5, -1376.4, -161.91};
	// Coefficients for c_2 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c2CoefficientsHigh = {
		834.91, 1981.8, 1885.7, 1027.1, 296.69, 29.902};
	// Coefficients for c_3 in the 2D E_f,HeV fit
	constexpr Kokkos::Array<double, 6> c3CoefficientsHigh = {
		1547.2, 3532.3, 3383.6, 1969.2, 695.17, 119.23};

	/**
	 * The formation energies for He_xV_1. The entry at i = 0 is for a single
	 * vacancy (He_0V_1) and is there as a buffer. Like the formation energies,
	 * i = heSize.
	 */
	constexpr Kokkos::Array<double, 15> heV1FormationEnergies = {0.0, 5.14166,
		8.20919, 11.5304, 14.8829, 18.6971, 22.2847, 26.3631, 30.1049, 34.0081,
		38.2069, 42.4217, 46.7378, 51.1551, 55.6738};

	/**
	 * The formation energies for He_xV_2. The entry at i = 0 is for a
	 * di-vacancy (He_0V_2) and is there as a buffer. Like the formation
	 * energies, i = heSize.
	 */
	constexpr Kokkos::Array<double, 19> heV2FormationEnergies = {0.0, 7.10098,
		8.39913, 9.41133, 11.8748, 14.8296, 17.7259, 20.7747, 23.7993, 26.7984,
		30.0626, 33.0385, 36.5173, 39.9406, 43.48, 46.8537, 50.4484, 54.0879,
		57.7939};

	// Initial declarations
	double energy = -util::infinity<double>;
	// The following coefficients are computed using the above and are used
	// to evaluate the full function f(x,y).
	Kokkos::Array<double, 4> coefficients{0.0, 0.0, 0.0, 0.0};

	// Check to see if the vacancy size is large enough that the energy can
	// be computed from the fit or if it is so small that the exact values
	// must be used instead.
	auto numHe = comp[Species::He];
	auto numV = comp[Species::V];
	if (numV > 2) {
		// Get the He/V ratio
		double x = 2.0 * (static_cast<double>(numHe) / (9.0 * numV)) - 1.0;
		// Initialize the vacancy number
		double y = 0.0;

		// We have 2 fits, one for low V and one for high V
		if (numV <= 27) {
			// Compute the vacancy number
			y = 2.0 * ((numV - 1.0) / 26.0) - 1.0;
			// Get the coefficients
			coefficients[0] =
				util::computeNthOrderLegendre<5>(x, c0CoefficientsLow);
			coefficients[1] =
				util::computeNthOrderLegendre<5>(x, c1CoefficientsLow);
			coefficients[2] =
				util::computeNthOrderLegendre<5>(x, c2CoefficientsLow);
			coefficients[3] =
				util::computeNthOrderLegendre<5>(x, c3CoefficientsLow);
		}
		else {
			// Compute the vacancy number
			y = 2.0 * ((numV - 1.0) / 451.0) - 1.0;
			// Get the coefficients
			coefficients[0] =
				util::computeNthOrderLegendre<5>(x, c0CoefficientsHigh);
			coefficients[1] =
				util::computeNthOrderLegendre<5>(x, c1CoefficientsHigh);
			coefficients[2] =
				util::computeNthOrderLegendre<5>(x, c2CoefficientsHigh);
			coefficients[3] =
				util::computeNthOrderLegendre<5>(x, c3CoefficientsHigh);
		}

		energy = util::computeNthOrderLegendre<3>(y, coefficients);
	}
	else if ((numV == 1 && numHe < heV1FormationEnergies.size()) ||
		(numV == 2 && numHe < heV2FormationEnergies.size())) {
		// Get the exact energy
		energy = (numV == 1) ? heV1FormationEnergies[numHe] :
							   heV2FormationEnergies[numHe];
	}

	// V Case
	if (numHe == 0 && numV < 3) {
		energy = (numV == 1) ? vFormationEnergies[0] : vFormationEnergies[1];
	}

	return energy;
}
} // namespace network
} // namespace core
} // namespace xolotl
