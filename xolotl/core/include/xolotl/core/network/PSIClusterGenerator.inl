#pragma once

#include <xolotl/util/MathUtils.h>
#include <xolotl/core/Constants.h>

namespace xolotl
{
namespace core
{
namespace network
{
PSIClusterGenerator<PSIFullSpeciesList>::PSIClusterGenerator(
        const options::IOptions& opts)
    :
    _hydrogenRadiusFactor(opts.getHydrogenFactor()),
    _maxHe(opts.getMaxImpurity()),
    _maxD(opts.getMaxD()),
    _maxT(opts.getMaxT()),
    _maxV(opts.getMaxV()),
    _groupingMin(opts.getGroupingMin()),
    _groupingWidthA(opts.getGroupingWidthA()),
    _groupingWidthB(opts.getGroupingWidthB())
{
}

PSIClusterGenerator<PSIFullSpeciesList>::PSIClusterGenerator(
        const options::IOptions& opts, std::size_t refineDepth)
    :
    Superclass(refineDepth),
    _hydrogenRadiusFactor(opts.getHydrogenFactor()),
    _maxHe(opts.getMaxImpurity()),
    _maxD(opts.getMaxD()),
    _maxT(opts.getMaxT()),
    _maxV(opts.getMaxV()),
    _groupingMin(opts.getGroupingMin()),
    _groupingWidthA(opts.getGroupingWidthA()),
    _groupingWidthB(opts.getGroupingWidthB())
{
}

KOKKOS_INLINE_FUNCTION
bool
PSIClusterGenerator<PSIFullSpeciesList>::intersect(const Region& region) const
{
    // I is never grouped
    if (region[Species::I].begin() > 0) {
        return true;
    }
    
    // V is never grouped
    if (region[Species::V].end() > 1 && region[Species::He].begin() == 0 &&
                region[Species::D].begin() == 0 && region[Species::T].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }
    
    //He is never grouped
    if (region[Species::He].end() > 1 && region[Species::D].begin() == 0 &&
                region[Species::T].begin() == 0 && region[Species::V].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }
    
    //D is never grouped
    if (region[Species::D].end() > 1 && region[Species::He].begin() == 0 &&
                region[Species::T].begin() == 0 && region[Species::V].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }
    
    //T is never grouped
    if (region[Species::T].end() > 1 && region[Species::He].begin() == 0 &&
                region[Species::D].begin() == 0 && region[Species::V].begin() == 0 &&
                region[Species::I].begin() == 0) {
        return true;
    }
    
    // Don't group under the given min for V
    if (region[Species::V].begin() < _groupingMin) {
        return true;
    }
    
    // Don't group above maxV so that the cluster are rejected by select
    if (region[Species::V].begin() >= _maxV) {
        return true;
    }
   
    // Else refine around the edge
    auto maxDPerV =
        KOKKOS_LAMBDA (AmountType amtV) { return (2.0/3.0) * getMaxHePerV(amtV) * (_maxD > 0); };
    auto maxTPerV =
        KOKKOS_LAMBDA (AmountType amtV) { return (2.0/3.0) * getMaxHePerV(amtV) * (_maxT > 0); };
    if (region[Species::V].end() > 1) {
        Composition lo = region.getOrigin();
        Composition hi = region.getUpperLimitPoint();
        
        if (lo[Species::He] <= getMaxHePerV(hi[Species::V]-1) 
        		&& hi[Species::He] - 1 >= getMaxHePerV(lo[Species::V] - 1)) {
        	return true;
        }
		if (lo[Species::D] <= maxDPerV(hi[Species::V]-1) 
        		&& hi[Species::D] - 1 >= maxDPerV(lo[Species::V] - 1) && _maxD > 0) {
            return true;
        }
		if (lo[Species::T] <= maxTPerV(hi[Species::V]-1) 
        		&& hi[Species::T] - 1 >= maxTPerV(lo[Species::V] - 1) && _maxT > 0) {
            return true;
        }
        auto hiH = hi[Species::D] + hi[Species::T];
        if (lo[Species::He] == 0) {
            return true;
        }
        if (hiH >= (2.0/3.0) * lo[Species::He] + 0.5 ) {
            return true;
        }
    }

    if (region[Species::V].length() == _groupingWidthB) {
        return false;
    }

    if (region[Species::He].length() == _groupingWidthA) {
        return false;
    }

    if (region[Species::D].length() == _groupingWidthA) {
        return false;
    }

    if (region[Species::T].length() == _groupingWidthA) {
        return false;
    }
    
    return true;
}

KOKKOS_INLINE_FUNCTION
bool
PSIClusterGenerator<PSIFullSpeciesList>::select(const Region& region) const
{
    // Remove 0
    if (region[Species::He].end() == 1 && region[Species::D].end() == 1 &&
            region[Species::T].end() == 1 && region[Species::V].end() == 1 &&
            region[Species::I].end() == 1) {
        return false;
    }

    // Interstitials
    if (region[Species::I].begin() > 0 && (region[Species::He].begin() > 0 ||
            region[Species::D].begin() > 0 || region[Species::T].begin() > 0 ||
            region[Species::V].begin() > 0)) {
        return false;
    }

    // Helium
    if (region[Species::He].begin() > _maxHe && region[Species::D].end() == 1 &&
            region[Species::T].end() == 1 && region[Species::V].end() == 1 &&
            region[Species::I].end() == 1) {
        return false;
    }

    // Deuterium
    if (region[Species::D].begin() > _maxD && region[Species::He].end() == 1 &&
            region[Species::T].end() == 1 && region[Species::V].end() == 1 &&
            region[Species::I].end() == 1) {
        return false;
    }

    // Tritium
    if (region[Species::T].begin() > _maxT && region[Species::He].end() == 1 &&
            region[Species::D].end() == 1 && region[Species::V].end() == 1 &&
            region[Species::I].end() == 1) {
        return false;
    }

    // Vacancy
    if (region[Species::V].begin() > _maxV) {
        return false;
    }

    // Can't cluster without V
    if (region[Species::V].end() == 1) {
        if (region[Species::He].begin() > 0 && region[Species::D].begin() > 0)
            return false;
        if (region[Species::He].begin() > 0 && region[Species::T].begin() > 0)
            return false;
        if (region[Species::D].begin() > 0 && region[Species::T].begin() > 0)
            return false;
    }

    // The edge
    auto maxDPerV =
        KOKKOS_LAMBDA (AmountType amtV) { return (2.0/3.0) * getMaxHePerV(amtV); };
    if (region[Species::V].end() > 1) {
        Composition lo = region.getOrigin();
        Composition hi = region.getUpperLimitPoint();

        //Too many helium
        if (lo[Species::He] > getMaxHePerV(hi[Species::V]-1)) {
            return false;
        }

        //Too many deuterium
        if (lo[Species::D] > maxDPerV(hi[Species::V]-1)) {
            return false;
        }

        //Too many tritium
        if (lo[Species::T] > maxDPerV(hi[Species::V]-1)) {
            return false;
        }

        // Too many hydrogen
        auto loH = lo[Species::D] + lo[Species::T];
        if (lo[Species::He] == 0) {
            if (loH > 6 * (hi[Species::V]-1)) {
                return false;
            }
        }
        else {
            if (loH > (2.0/3.0) * (hi[Species::He] - 1) + 0.5 ) {
                return false;
            }
        }
    }

    return true;
}

KOKKOS_INLINE_FUNCTION
typename PSIClusterGenerator<PSIFullSpeciesList>::AmountType
PSIClusterGenerator<PSIFullSpeciesList>::getMaxHePerV(AmountType amtV) noexcept
{
    /**
     * The maximum number of helium atoms that can be combined with a
     * vacancy cluster with size equal to the index i in the array plus one.
     * For example, an HeV size cluster with size 1 would have
     * size = i+1 = 1 and i = 0. It could support a mixture of up to nine
     * helium atoms with one vacancy.
     */
    constexpr Kokkos::Array<AmountType, 30> maxHePerV = {
        0, 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
        85, 90, 95, 98, 100, 101, 103, 105, 107, 109, 110, 112, 116
    };

    if (amtV < maxHePerV.size()) {
        return maxHePerV[amtV];
    }
    return 4 * amtV;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<PSIFullSpeciesList>::getFormationEnergy(
    const Cluster<PlsmContext>& cluster) const noexcept
{
    constexpr auto infinity = std::numeric_limits<double>::infinity();

    // I formation energies in eV
    constexpr Kokkos::Array<double, 7> iFormationEnergies = {
        0.0, 10.0, 18.5, 27.0, 35.0, 42.5, 48.0
    };
    // He formation energies in eV
    constexpr Kokkos::Array<double, 9> heFormationEnergies = {
        0.0, 6.15, 11.44, 16.35, 21.0, 26.1, 30.24, 34.93, 38.80
    };


    const auto& reg = cluster.getRegion();
    double formationEnergy {};
    // All the possibly grouped ones are 0
    if (reg[Species::V].end() >= _groupingMin) return formationEnergy;
    
    if (reg.isSimplex()) {
        Composition comp(reg.getOrigin());
        if (comp.isOnAxis(Species::I)) {
            auto amtI = comp[Species::I];
            if (amtI < iFormationEnergies.size()) {
                formationEnergy = iFormationEnergies[amtI];
            }
            else {
                formationEnergy = /* 48 + 6*(amtI - 6) */
                    6.0 * (2.0 + amtI);
            }
        }
        else if (comp.isOnAxis(Species::He)) {
            auto amtHe = comp[Species::He];
            if (amtHe < heFormationEnergies.size()) {
                formationEnergy = heFormationEnergies[amtHe];
            }
            else {
                formationEnergy = infinity;
            }
        }
        else if (comp.isOnAxis(Species::D)) {
            formationEnergy = infinity;
        }
        else if (comp.isOnAxis(Species::T)) {
            formationEnergy = infinity;
        }
        else {
            formationEnergy = getHeVFormationEnergy(comp);
        }
    }
    return formationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<PSIFullSpeciesList>::getMigrationEnergy(
    const Cluster<PlsmContext>& cluster) const noexcept
{
    // I migration energies in eV
    constexpr Kokkos::Array<double, 6> iMigration = {
        0.0, 0.01, 0.02, 0.03, 0.04, 0.05
    };
    // He migration energies in eV
    constexpr Kokkos::Array<double, 8> heMigration = {
        0.0, 0.13, 0.20, 0.25, 0.20, 0.12, 0.3, 0.4
    };
    // The migration energy for a single deuterium.
    constexpr double dOneMigrationEnergy = 0.38;
    // The migration energy for a single tritium.
    static constexpr double tOneMigrationEnergy = 0.38;
    // The migration energy for a single vacancy in eV
    static constexpr double vOneMigration = 1.30;

    const auto& reg = cluster.getRegion();
    double migrationEnergy = std::numeric_limits<double>::infinity();
    if (reg.isSimplex()) {
        Composition comp(reg.getOrigin());
        if (comp.isOnAxis(Species::I)) {
            auto amtI = comp[Species::I];
            if (amtI < iMigration.size()) {
                migrationEnergy = iMigration[amtI];
            }
            else {
                migrationEnergy = util::min((double) amtI, 15.0) * 0.1;
            }
        }
        else if (comp.isOnAxis(Species::He)) {
            auto amtHe = comp[Species::He];
            if (amtHe < heMigration.size()) {
                migrationEnergy = heMigration[amtHe];
            }
        }
        else if (comp.isOnAxis(Species::D)) {
            if (comp[Species::D] == 1) {
                migrationEnergy = dOneMigrationEnergy;
            }
        }
        else if (comp.isOnAxis(Species::T)) {
            if (comp[Species::T] == 1) {
                migrationEnergy = tOneMigrationEnergy;
            }
        }
        else if (comp.isOnAxis(Species::V)) {
            if (comp[Species::V] <= 1) {
                migrationEnergy = vOneMigration;
            }
        }
    }
    return migrationEnergy;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<PSIFullSpeciesList>::getDiffusionFactor(
    const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
    // I diffusion factors in nm^2/s
    constexpr Kokkos::Array<double, 6> iDiffusion = {
        0.0, 8.8e+10, 8.0e+10, 3.9e+10, 2.0e+10, 1.0e+10
    };
    // He diffusion factors in nm^2/s
    constexpr Kokkos::Array<double, 8> heDiffusion = {
        0.0, 2.9e+10, 3.2e+10, 2.3e+10, 1.7e+10, 5.0e+09, 1.0e+09, 5.0e+08
    };
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
                diffusionFactor = iDiffusion[amtI];
            }
            else {
                diffusionFactor = iDiffusion[1] / (double) amtI;
            }
        }
        else if (comp.isOnAxis(Species::He)) {
            auto amtHe = comp[Species::He];
            if (amtHe < heDiffusion.size()) {
                diffusionFactor = heDiffusion[amtHe];
            }
        }
        else if (comp.isOnAxis(Species::D)) {
            if (comp[Species::D] == 1) {
                diffusionFactor = dOneDiffusionFactor;
            }
        }
        else if (comp.isOnAxis(Species::T)) {
            if (comp[Species::T] == 1) {
                diffusionFactor = tOneDiffusionFactor;
            }
        }
        else if (comp.isOnAxis(Species::V)) {
            if (comp[Species::V] == 1) {
                diffusionFactor = vOneDiffusion;
            }
        }
    }

    return diffusionFactor;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<PSIFullSpeciesList>::getReactionRadius(
    const Cluster<PlsmContext>& cluster,
    double latticeParameter, double interstitialBias, double impurityRadius)
    const noexcept
{
    const auto& reg = cluster.getRegion();
    double radius = 0.0;
    if (reg.isSimplex()) {
        Composition comp(reg.getOrigin());
        if (comp.isOnAxis(Species::I)) {
            double EightPi = 8.0 * ::xolotl::core::pi;
            double aCubed = pow(latticeParameter, 3.0);
            double termOne = interstitialBias * (sqrt(3.0) / 4.0)
                    * latticeParameter;
            double termTwo = pow((3.0 / EightPi) * aCubed * comp[Species::I], (1.0 / 3.0));
            double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));
            radius = termOne + termTwo - termThree;
        }
        else if (comp.isOnAxis(Species::He)) {
            double FourPi = 4.0 * ::xolotl::core::pi;
            double aCubed = pow(latticeParameter, 3);
            double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::He],
                    (1.0 / 3.0));
            double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed,
                    (1.0 / 3.0));
            radius = impurityRadius + termOne - termTwo;
        }
        else if (comp.isOnAxis(Species::D)) {
            double FourPi = 4.0 * ::xolotl::core::pi;
            double aCubed = pow(latticeParameter, 3);
            double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::D],
                    (1.0 / 3.0));
            double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed,
                    (1.0 / 3.0));
            radius = (impurityRadius + termOne - termTwo) * _hydrogenRadiusFactor;
        }
        else if (comp.isOnAxis(Species::T)) {
            double FourPi = 4.0 * ::xolotl::core::pi;
            double aCubed = pow(latticeParameter, 3);
            double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * comp[Species::T],
                    (1.0 / 3.0));
            double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed,
                    (1.0 / 3.0));
            radius = (impurityRadius + termOne - termTwo) * _hydrogenRadiusFactor;
        }
        else {
            radius = (sqrt(3.0) / 4.0) * latticeParameter
                    + pow(
                            (3.0 * pow(latticeParameter, 3.0) * comp[Species::V])
                                    / (8.0 * ::xolotl::core::pi), (1.0 / 3.0))
                    - pow((3.0 * pow(latticeParameter, 3.0)) / (8.0 * ::xolotl::core::pi),
                            (1.0 / 3.0));
        }
    }
    else {
        // Loop on the V range
        for (auto j : makeIntervalRange(reg[Species::V])) {
            radius += (sqrt(3.0) / 4.0) * latticeParameter
                    + pow(
                            (3.0 * pow(latticeParameter, 3.0) * (double) j)
                                    / (8.0 * ::xolotl::core::pi), (1.0 / 3.0))
                    - pow((3.0 * pow(latticeParameter, 3.0)) / (8.0 * ::xolotl::core::pi),
                            (1.0 / 3.0));
        }
        // Average the radius
        radius /= reg[Species::V].length();
    }

    return radius;
}

KOKKOS_INLINE_FUNCTION
double
PSIClusterGenerator<PSIFullSpeciesList>::getHeVFormationEnergy(Composition comp)
    const noexcept
{
    // V formation energies in eV
    constexpr Kokkos::Array<double, 2> vFormationEnergies = { 3.6, 7.25 };

    // Coefficients for the Legendre polynomial fit
    // Low means V <= 27
    // Coefficients for c_0 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c0CoefficientsLow = {
        253.35, 435.36, 336.50, 198.92, 95.154, 21.544
    };
    // Coefficients for c_1 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c1CoefficientsLow = {
        493.29, 1061.3, 1023.9, 662.92, 294.24, 66.962
    };
    // Coefficients for c_2 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c2CoefficientsLow = {
        410.40, 994.89, 1044.6, 689.41, 286.52, 60.712
    };
    // Coefficients for c_3 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c3CoefficientsLow = {
        152.99, 353.16, 356.10, 225.75, 87.077, 15.640
    };
    // High means V > 27
    // Coefficients for c_0 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c0CoefficientsHigh = {
        -847.90, -3346.9, -4510.3, -3094.7, -971.18, -83.770
    };
    // Coefficients for c_1 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c1CoefficientsHigh = {
        -1589.3, -4894.6, -6001.8, -4057.5, -1376.4, -161.91
    };
    // Coefficients for c_2 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c2CoefficientsHigh = {
        834.91, 1981.8, 1885.7, 1027.1, 296.69, 29.902
    };
    // Coefficients for c_3 in the 2D E_f,HeV fit
    constexpr Kokkos::Array<double, 6> c3CoefficientsHigh = {
        1547.2, 3532.3, 3383.6, 1969.2, 695.17, 119.23
    };

    /**
     * The formation energies for He_xV_1. The entry at i = 0 is for a single
     * vacancy (He_0V_1) and is there as a buffer. Like the formation energies,
     * i = heSize.
     */
    constexpr Kokkos::Array<double, 15> heV1FormationEnergies = {
        0.0, 5.14166, 8.20919, 11.5304, 14.8829, 18.6971, 22.2847, 26.3631,
        30.1049, 34.0081, 38.2069, 42.4217, 46.7378, 51.1551, 55.6738
    };

    /**
     * The formation energies for He_xV_2. The entry at i = 0 is for a
     * di-vacancy (He_0V_2) and is there as a buffer. Like the formation
     * energies, i = heSize.
     */
    constexpr Kokkos::Array<double, 19> heV2FormationEnergies = {
        0.0, 7.10098, 8.39913, 9.41133, 11.8748, 14.8296, 17.7259, 20.7747,
        23.7993, 26.7984, 30.0626, 33.0385, 36.5173, 39.9406, 43.48,
        46.8537, 50.4484, 54.0879, 57.7939
    };

    // Initial declarations
    double energy = -std::numeric_limits<double>::infinity();
    // The following coefficients are computed using the above and are used
    // to evaluate the full function f(x,y).
    Kokkos::Array<double, 4> coefficients { 0.0, 0.0, 0.0, 0.0 };

    // Check to see if the vacancy size is large enough that the energy can
    // be computed from the fit or if it is so small that the exact values
    // must be used instead.
    auto numHe = comp[Species::He];
    auto numV = comp[Species::V];
    if (numV > 2) {
        // Get the He/V ratio
        double x = 2.0 * (static_cast<double>(numHe) / (9.0*numV)) - 1.0;
        // Initialize the vacancy number
        double y = 0.0;

        // We have 2 fits, one for low V and one for high V
        if (numV <= 27) {
            // Compute the vacancy number
            y = 2.0 * ((numV - 1.0) / 26.0) - 1.0;
            // Get the coefficients
            coefficients[0] = util::computeNthOrderLegendre<5>(x, c0CoefficientsLow);
            coefficients[1] = util::computeNthOrderLegendre<5>(x, c1CoefficientsLow);
            coefficients[2] = util::computeNthOrderLegendre<5>(x, c2CoefficientsLow);
            coefficients[3] = util::computeNthOrderLegendre<5>(x, c3CoefficientsLow);
        }
        else {
            // Compute the vacancy number
            y = 2.0 * ((numV - 1.0) / 451.0) - 1.0;
            // Get the coefficients
            coefficients[0] = util::computeNthOrderLegendre<5>(x, c0CoefficientsHigh);
            coefficients[1] = util::computeNthOrderLegendre<5>(x, c1CoefficientsHigh);
            coefficients[2] = util::computeNthOrderLegendre<5>(x, c2CoefficientsHigh);
            coefficients[3] = util::computeNthOrderLegendre<5>(x, c3CoefficientsHigh);
        }

        energy = util::computeNthOrderLegendre<3>(y, coefficients);
    }
    else if ((numV == 1 && numHe < heV1FormationEnergies.size()) ||
            (numV == 2 && numHe < heV2FormationEnergies.size())) {
        // Get the exact energy
        energy = (numV == 1) ?
            heV1FormationEnergies[numHe] :
            heV2FormationEnergies[numHe];
    }

    // V Case
    if (numHe == 0 && numV < 3) {
        energy = (numV == 1) ? vFormationEnergies[0] : vFormationEnergies[1];
    }

    return energy;
}
}
}
}
