#pragma once

#include <MathUtils.h>

namespace xolotlCore
{
namespace experimental
{
KOKKOS_INLINE_FUNCTION
bool
FeClusterGenerator<FeFullSpeciesList>::intersect(const Region& region) const
{
    // I is never grouped
    if (region[Species::I].begin() > 0) {
        return true;
    }

    // V is never grouped
    if (region[Species::V].end() > 1 && region[Species::He].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }

    //He is never grouped
    if (region[Species::He].end() > 1 && region[Species::V].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }

    return true;
}

KOKKOS_INLINE_FUNCTION
bool
FeClusterGenerator<FeFullSpeciesList>::select(const Region& region) const
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

    // Vacancy
    if (region[Species::V].begin() > 10 && region[Species::He].end() == 1 &&
            region[Species::I].end() == 1) {
        return false;
    }

    return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
FeClusterGenerator<FeFullSpeciesList>::getMigrationEnergy(
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
FeClusterGenerator<FeFullSpeciesList>::getDiffusionFactor(
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
FeClusterGenerator<FeFullSpeciesList>::getReactionRadius(
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
