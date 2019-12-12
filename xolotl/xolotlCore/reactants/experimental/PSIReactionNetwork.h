#pragma once

#include <MathUtils.h>

#include <experimental/ReactionNetwork.h>
#include <experimental/PSIReaction.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;


template <typename TSpeciesEnum>
class PSIReaction;


template <typename TSpeciesEnum>
class PSIClusterGenerator
{
};


enum class PSIFullSpeciesList
{
    He,
    D,
    T,
    V,
    I
};


template <>
struct hasInterstitial<PSIFullSpeciesList> : std::true_type { };


template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 5;

    using ReactionType = PSIReaction<Species>;

    using ClusterGenerator = PSIClusterGenerator<Species>;
};


template <typename TSpeciesEnum>
class PSIReactionNetwork :
    public ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>
{
    friend class ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
public:
    using Superclass = ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
    using Subpaving = typename Superclass::Subpaving;
    using ReactionType = typename Superclass::ReactionType;
    using Composition = typename Superclass::Composition;
    using Species = PSIFullSpeciesList;

    using Superclass::Superclass;

private:
    double
    checkLatticeParameter(double latticeParameter)
    {
        if (latticeParameter <= 0.0) {
            return tungstenLatticeConstant;
        }
        return latticeParameter;
    }

    double
    checkImpurityRadius(double impurityRadius)
    {
        if (impurityRadius <= 0.0) {
            return heliumRadius;
        }
        return impurityRadius;
    }

    typename Superclass::ClusterSetsPair
    defineReactionClusterSets(
        typename Subpaving::template TilesView<plsm::OnHost> tiles,
        Kokkos::View<double*>::HostMirror diffusionFactor)
    {
        typename Superclass::ClusterSetsPair ret;

        auto species = this->getSpeciesRange();
        auto speciesNoI = this->getSpeciesRangeNoI();

        std::size_t numClusters = tiles.extent(0);
        for (std::size_t i = 0; i < numClusters; ++i)
            for (std::size_t j = i; j < numClusters; ++j) {
                // Check that at least one of them is mobile
                if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0)
                    continue;

                // Get the composition of each cluster
                const auto& cl1Reg = tiles(i).getRegion();
                const auto& cl2Reg = tiles(j).getRegion();
                Composition lo1 = cl1Reg.getOrigin();
                Composition hi1 = cl1Reg.getUpperLimitPoint();
                Composition lo2 = cl2Reg.getOrigin();
                Composition hi2 = cl2Reg.getUpperLimitPoint();

                // Special case for I + I
                if (cl1Reg.isSimplex() && cl2Reg.isSimplex() && lo1.isOnAxis(Species::I)
                    && lo2.isOnAxis(Species::I)) {
                    // Compute the composition of the new cluster
                    std::size_t size = lo1[Species::I]
                        + lo2[Species::I];
                    // Loop on potential products
                    for (std::size_t k = 0; k < numClusters; ++i) {
                        // Get the composition
                        const auto& prodReg = tiles(k).getRegion();
                        if (!prodReg.isSimplex()) continue;
                        Composition lo = prodReg.getOrigin();
                        if (!lo.isOnAxis(Species::I)) continue;
                        if (lo[Species::I] == size) {
                            ret.prodClusterSets.emplace_back(i, j, k);
                            if (lo1[Species::I] == 1
                                || lo2[Species::I] == 1)
                                ret.dissClusterSets.emplace_back(k, i, j);
                            break;
                        }
                    }
                    continue;
                }

                // Special case for I + V
                if (cl1Reg.isSimplex() && cl2Reg.isSimplex() && ((lo1.isOnAxis(Species::I)
                        && lo2.isOnAxis(Species::V)) || (lo1.isOnAxis(Species::V)
                                && lo2.isOnAxis(Species::I)))) {
                    // Find out which one is which
                    auto vSize = lo1.isOnAxis(Species::V)? lo1[Species::V] :
                       lo2[Species::V];
                    auto iSize = lo1.isOnAxis(Species::I)? lo1[Species::I] :
                       lo2[Species::I];
                    // Compute the product size
                    int prodSize = vSize - iSize;
                    // 3 cases
                    if (prodSize > 0) {
                        // Looking for V cluster
                        for (std::size_t k = 0; k < numClusters; ++i) {
                            // Get the composition
                            const auto& prodReg = tiles(k).getRegion();
                            if (!prodReg.isSimplex()) continue;
                            Composition lo = prodReg.getOrigin();
                            if (!lo.isOnAxis(Species::V)) continue;
                            if (lo[Species::V] == prodSize) {
                                ret.prodClusterSets.emplace_back(i, j, k);
                                // No dissociation
                                break;
                            }
                        }
                    }
                    else if (prodSize < 0) {
                        // Looking for I cluster
                        for (std::size_t k = 0; k < numClusters; ++i) {
                            // Get the composition
                            const auto& prodReg = tiles(k).getRegion();
                            if (!prodReg.isSimplex()) continue;
                            Composition lo = prodReg.getOrigin();
                            if (!lo.isOnAxis(Species::I)) continue;
                            if (lo[Species::I] == prodSize) {
                                ret.prodClusterSets.emplace_back(i, j, k);
                                // No dissociation
                                break;
                            }
                        }
                    }
                    else {
                        // No product
                        ret.prodClusterSets.emplace_back(i, j);
                    }
                    continue;
                }

                // General case
                std::vector<std::pair<int, int> > bounds;
                // Loop on the species
                for (auto l : species) {
                    auto low = lo1[l] + lo2[l];
                    auto high = hi1[l] + hi2[l] - 2;
                    // Special case for I
                    if (l == Species::I) {
                        // Look for V
                        size_t vIdx = 0;
                        for (auto m : speciesNoI) {
                            if (m == Species::V) {
                                vIdx = m();
                                break;
                            }
                        }
                        bounds[vIdx].first -= high;
                        bounds[vIdx].second -= low;

                        if (bounds[vIdx].first < 0) bounds[vIdx].first = 0;
                        if (bounds[vIdx].second < 0) bounds[vIdx].second = 0;
                    }
                    else
                        bounds.emplace_back(std::make_pair(low, high));
                }

                // Look for potential product
                std::size_t nProd = 0;
                for (std::size_t k = 0; k < numClusters; ++i) {
                    // Get the composition
                    const auto& prodReg = tiles(k).getRegion();
                    bool isGood = true;
                    // Loop on the species
                    // TODO: check l correspond to the same species in bounds and prod
                    for (auto l : speciesNoI) {
                        if (prodReg[l()].begin() > bounds[l()].second) {
                            isGood = false;
                            break;
                        }
                        if (prodReg[l()].end() - 1 < bounds[l()].first) {
                            isGood = false;
                            break;
                        }
                    }

                    if (isGood) {
                        // Increase nProd
                        nProd++;
                        ret.prodClusterSets.emplace_back(i, j, k);
                        // TODO: will have to add some rules, i or j should be a simplex cluster of max size 1
                        ret.dissClusterSets.emplace_back(k, i, j);
                    }
                }

                // Special case for trap-mutation
                if (nProd == 0) {
                    // Look for larger clusters
                    // Loop on possible I sizes
                    // TODO: get the correct value for maxISize
                    std::size_t maxISize = 6;
                    for (std::size_t n = 1; n <= maxISize; ++n) {
                        // Find the corresponding cluster
                        Composition comp{};
                        comp[Species::I] = n;
                        auto iCluster = this->findCluster(comp);

                        // Look for V
                        size_t vIdx = 0;
                        for (auto m : speciesNoI) {
                            if (m == Species::V) {
                                vIdx = m();
                                break;
                            }
                        }
                        bounds[vIdx].first += 1;
                        bounds[vIdx].second += 1;

                        // Look for potential product
                        std::size_t nProd = 0;
                        for (std::size_t k = 0; k < numClusters; ++i) {
                            // Get the composition
                            const auto& prodReg = tiles(k).getRegion();
                            bool isGood = true;
                            // Loop on the species
                            // TODO: check l correspond to the same species in bounds and prod
                            for (auto l : speciesNoI) {
                                if (prodReg[l()].begin() > bounds[l()].second) {
                                    isGood = false;
                                    break;
                                }
                                if (prodReg[l()].end() - 1 < bounds[l()].first) {
                                    isGood = false;
                                    break;
                                }
                            }

                            if (isGood) {
                                // Increase nProd
                                nProd++;
                                ret.prodClusterSets.emplace_back(i, j, k, iCluster.getId());
                                // No dissociation
                            }
                        }
                        // Stop if we found a product
                        if (nProd > 0) break;
                    }
                }
        }

        return ret;
    }
};


template <>
class PSIClusterGenerator<PSIFullSpeciesList> :
    public
    plsm::refine::Detector<PSIClusterGenerator<PSIFullSpeciesList>>
{
public:
    using Species = PSIFullSpeciesList;
    using Superclass = plsm::refine::Detector<PSIClusterGenerator<Species>>;
    using NetworkType = PSIReactionNetwork<Species>;
    using Cluster = typename NetworkType::Cluster;
    using Region = typename NetworkType::Region;
    using Composition = typename NetworkType::Composition;
    using AmountType = typename NetworkType::AmountType;

    PSIClusterGenerator(const IOptions& options)
        :
        _hydrogenRadiusFactor(options.getHydrogenFactor())
    {
    }

    PSIClusterGenerator(const IOptions& options, std::size_t refineDepth)
        :
        Superclass(refineDepth),
        _hydrogenRadiusFactor(options.getHydrogenFactor())
    {
    }

    KOKKOS_INLINE_FUNCTION
    bool
    intersect(const Region& region) const
    {
        //TODO
        return true;
    }

    KOKKOS_INLINE_FUNCTION
    bool
    select(const Region& region) const
    {
        auto maxDPerV =
            KOKKOS_LAMBDA (AmountType amtV) { return (8.0/3.0) * amtV; };
        if (region[Species::V].begin() > 0) {
            Composition lo = region.getOrigin();
            Composition hi = region.getUpperLimitPoint();

            //Too many helium
            if (lo[Species::He] > getMaxHePerV(lo[Species::V]) &&
                    lo[Species::He] > getMaxHePerV(hi[Species::V] - 1)) {
                return false;
            }

            //Too many deuterium
            if (lo[Species::D] > maxDPerV(lo[Species::V]) &&
                    lo[Species::D] > maxDPerV(hi[Species::V])) {
                return false;
            }

            //Too many tritium
            if (lo[Species::T] > maxDPerV(lo[Species::V]) &&
                    lo[Species::T] > maxDPerV(hi[Species::V])) {
                return false;
            }

            // Too many hydrogen
            auto loH = lo[Species::D] + lo[Species::T];
            if (lo[Species::He] == 0 && hi[Species::He] - 1 == 0) {
                if (loH > 6 * lo[Species::V] && loH > 6 * hi[Species::V]) {
                    return false;
                }
            }
            else {
                if (loH > (2.0/3.0) * lo[Species::He] + 0.5 &&
                        loH > (2.0/3.0) * hi[Species::He] + 0.5) {
                    return false;
                }
            }
        }

        return true;
    }

    KOKKOS_INLINE_FUNCTION
    AmountType
    getMaxHePerV(AmountType amtV) const noexcept
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

    double
    getFormationEnergy(const Cluster& cluster) const noexcept
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

    double
    getMigrationEnergy(const Cluster& cluster) const noexcept
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

    double
    getDiffusionFactor(const Cluster& cluster) const noexcept
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
                if (comp[Species::V] <= 1) {
                    diffusionFactor = vOneDiffusion;
                }
            }
        }
        return diffusionFactor;
    }

private:
    //TODO: Move to MathUtils.h
    template <std::size_t N>
    static double
    computeNthOrderLegendre(double x, const Kokkos::Array<double, N+1>& coeffs)
    {
        int currDegree = 0;
        double value = 0.0;
        for (auto currCoeff : coeffs) {
            value += currCoeff * legendrePolynomial(x, currDegree);
            ++currDegree;
        }
        return value;
    }

    double
    getHeVFormationEnergy(Composition comp) const noexcept
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
                coefficients[0] = computeNthOrderLegendre<5>(x, c0CoefficientsLow);
                coefficients[1] = computeNthOrderLegendre<5>(x, c1CoefficientsLow);
                coefficients[2] = computeNthOrderLegendre<5>(x, c2CoefficientsLow);
                coefficients[3] = computeNthOrderLegendre<5>(x, c3CoefficientsLow);
            }
            else {
                // Compute the vacancy number
                y = 2.0 * ((numV - 1.0) / 451.0) - 1.0;
                // Get the coefficients
                coefficients[0] = computeNthOrderLegendre<5>(x, c0CoefficientsHigh);
                coefficients[1] = computeNthOrderLegendre<5>(x, c1CoefficientsHigh);
                coefficients[2] = computeNthOrderLegendre<5>(x, c2CoefficientsHigh);
                coefficients[3] = computeNthOrderLegendre<5>(x, c3CoefficientsHigh);
            }

            energy = computeNthOrderLegendre<3>(y, coefficients);
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

private:
    // The factor between He and H radius sizes
    double _hydrogenRadiusFactor {0.25};
};
}
}
