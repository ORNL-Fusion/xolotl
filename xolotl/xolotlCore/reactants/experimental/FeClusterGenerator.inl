#pragma once

#include <MathUtils.h>

namespace xolotlCore
{
namespace experimental
{
KOKKOS_INLINE_FUNCTION
bool
FeClusterGenerator::intersect(const Region& region) const
{
    // I is never grouped
    if (region[Species::I].begin() > 0) {
        return true;
    }

    // V is never grouped
    if (region[Species::V].end() > 1 && region[Species::V].begin() < 11 && region[Species::He].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }

    //He is never grouped
    if (region[Species::He].end() > 1 && region[Species::He].begin() < 9 && region[Species::V].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }
    
    // HeV
    if (region[Species::He].begin() < _groupingMin && region[Species::V].begin() < _groupingMin) {
        return true;
    }
    Composition lo = region.getOrigin();
    Composition hi = region.getUpperLimitPoint();
    double amtHe = (double) (lo[Species::He] + hi[Species::He] - 1) / 2.0,
        amtV = (double) (lo[Species::V] + hi[Species::V] - 1) / 2.0;
    double amt = sqrt(pow(amtHe, 2.0) + pow(amtV, 2.0));
    double ibe = 4.88 + 2.59
            * (pow(amtV, 2.0 / 3.0)
            - pow(amtV - 1.0, 2.0 / 3.0))
            - 2.5 * log(1.0 + (amtHe / amtV));
    auto distance = abs(ibe - 1.0);
    if (distance * 0.2 < 1.0) {
    if (region[Species::He].begin() < _groupingMin) {
        if (region[Species::V].length() < max(_groupingWidthV + 1, (AmountType) ((region[Species::V].begin() - _groupingMin) * 0.1))) {
            return false;
        }
        else return true;
    }
    if (region[Species::V].begin() < _groupingMin) {
        if (region[Species::He].length() < max(_groupingWidthHe + 1, (AmountType) ((region[Species::He].begin() - _groupingMin) * 0.1))) {
            return false;
        }
        else return true;
    }
    if (region[Species::He].length() < max(_groupingWidthHe + 1, (AmountType) ((amt - _groupingMin) * 0.1)) 
            || region[Species::V].length() < max(_groupingWidthV + 1, (AmountType) ((amt - _groupingMin) * 0.1))) {
        return false;
    }
    }
    else {
    if (region[Species::He].length() < max(_groupingWidthHe + 1, (AmountType) exp(distance * 1.0) * _groupingWidthHe * 2)
            || region[Species::V].length() < max(_groupingWidthV + 1, (AmountType) exp(distance * 1.0) * _groupingWidthV * 2)) return false;
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
    if (region[Species::I].begin() > 0 && (region[Species::He].begin() > 0 ||
            region[Species::V].begin() > 0)) {
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
    constexpr Kokkos::Array<double, 4> heMigration = {
        0.0, 0.06, 0.06, 0.06
    };
    // V migration energies in eV
    constexpr Kokkos::Array<double, 5> vMigration = {
        0.0, 0.67, 0.62, 0.37, 0.48
    };

    const auto& reg = cluster.getRegion();
    double migrationEnergy = std::numeric_limits<double>::infinity();
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
    const Cluster<PlsmContext>& cluster) const noexcept
{
    // I diffusion factors in nm^2/s
    constexpr double iOneDiffusionFactor = 1.0e+11;
    // He diffusion factors in nm^2/s
    constexpr Kokkos::Array<double, 4> heDiffusion = {
        0.0, 1.0e+11, 5.0e+10, 3.3e+10
    };
    // V diffusion factors in nm^2/s
    constexpr Kokkos::Array<double, 5> vDiffusion = {
        0.0, 1.0e+11, 5.0e+10, 3.3e+10, 2.5e+10
    };

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
FeClusterGenerator::getReactionRadius(
    const Cluster<PlsmContext>& cluster,
    double latticeParameter, double interstitialBias, double impurityRadius)
    const noexcept
{
    const auto& reg = cluster.getRegion();
    double radius = 0.0;
    if (reg.isSimplex()) {
        Composition comp(reg.getOrigin());
        if (comp.isOnAxis(Species::I)) {
            radius = latticeParameter * cbrt(3.0 / xolotlCore::pi) * 0.5;
        }
        else if (comp.isOnAxis(Species::He)) {
            double FourPi = 4.0 * xolotlCore::pi;
            double aCubed = pow(latticeParameter, 3);
            double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::He],
                    (1.0 / 3.0));
            double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed,
                    (1.0 / 3.0));
            radius = impurityRadius + termOne - termTwo;
        }
        else {
            radius = latticeParameter
                    * pow((3.0 * comp[Species::V]) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;
        }
    }
    else {
        // Loop on the V range
        for (auto j : makeIntervalRange(reg[Species::V])) {
            radius += latticeParameter
                            * pow((3.0 * (double) j) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;
        }
        // Average the radius
        radius /= reg[Species::V].length();
    }

    return radius;
}
}
}
