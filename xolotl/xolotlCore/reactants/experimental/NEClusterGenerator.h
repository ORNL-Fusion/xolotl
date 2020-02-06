#pragma once

#include <plsm/refine/Detector.h>

namespace xolotlCore
{
namespace experimental
{
class NEClusterGenerator : public plsm::refine::Detector<NEClusterGenerator>
{
public:
    using Superclass = plsm::refine::Detector<NEClusterGenerator>;
    using NetworkType = NEReactionNetwork;
    using Species = NESpecies;

    template <typename PlsmContext>
    using Cluster = Cluster<NetworkType, PlsmContext>;

    using Region = typename NetworkType::Region;
    using Composition = typename NetworkType::Composition;
    using AmountType = typename NetworkType::AmountType;

    NEClusterGenerator(const IOptions& options)
        :
        _maxXe(options.getMaxImpurity()),
        _xeDiffusivity(options.getXenonDiffusivity()),
        _xeDiffusive(_xeDiffusivity > 0.0),
        _groupingMin(options.getGroupingMin()),
        _groupingWidth(options.getGroupingWidthA()),
        _density(options.getDensity())
    {
    }

    NEClusterGenerator(const IOptions& options, std::size_t refineDepth)
        :
        Superclass(refineDepth),
        _maxXe(options.getMaxImpurity()),
        _xeDiffusivity(options.getXenonDiffusivity()),
        _xeDiffusive(_xeDiffusivity > 0.0),
        _groupingMin(options.getGroupingMin()),
        _groupingWidth(options.getGroupingWidthA()),
        _density(options.getDensity())
    {
    }

    KOKKOS_INLINE_FUNCTION
    bool
    intersect(const Region& region) const
    {
        if (region[Species::Xe].end() < _groupingMin) {
            return true;
        }
        if (region[Species::Xe].length() == _groupingWidth) {
            return false;
        }
        return true;
    }

    KOKKOS_INLINE_FUNCTION
    bool
    select(const Region& region) const
    {
        // Remove 0
        if (region[Species::Xe].end() == 1) {
            return false;
        }

        return true;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getFormationEnergy(const Cluster<PlsmContext>& cluster) const noexcept
    {
        /**
         * The set of xenon formation energies up to Xe_29 indexed by size. That is
         * E_(f,Xe_1) = xeFormationEnergies[1]. The value at index zero is just
         * padding to make the indexing easy.
         */
        constexpr Kokkos::Array<double, 30> xeFormationEnergies = {
            0.0, 7.0, 12.15, 17.15, 21.90, 26.50, 31.05, 35.30, 39.45, 43.00,
            46.90, 50.65, 53.90, 56.90, 59.80, 62.55, 65.05, 67.45, 69.45,
            71.20, 72.75, 74.15, 75.35, 76.40, 77.25, 77.95, 78.45, 78.80,
            78.95, 79.0
        };

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
        else return 79.0;
        return 0.0;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getMigrationEnergy(const Cluster<PlsmContext>& cluster) const noexcept
    {
        constexpr double xeOneMigration = 0.0;

        const auto& reg = cluster.getRegion();
        double migrationEnergy = std::numeric_limits<double>::infinity();
        if (reg.isSimplex()) {
            auto amtXe = reg.getOrigin()[0];
            if (amtXe <= 1) {
                migrationEnergy = _xeDiffusive ? 0.0 : xeOneMigration;
            }
        }
        return migrationEnergy;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getDiffusionFactor(const Cluster<PlsmContext>& cluster) const noexcept
    {
        constexpr double xeOneDiffusion = 1.0;

        const auto& reg = cluster.getRegion();
        double diffusionFactor = 0.0;
        if (reg.isSimplex()) {
            auto amtXe = reg.getOrigin()[0];
            if (amtXe == 1) {
                diffusionFactor = _xeDiffusive ? _xeDiffusivity : xeOneDiffusion;
            }
        }
        return diffusionFactor;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getReactionRadius(const Cluster<PlsmContext>& cluster,
        double latticeParameter, double interstitialBias, double impurityRadius)
        const noexcept
    {
        const auto& reg = cluster.getRegion();
        double radius = 0.0;
        double FourPi = 4.0 * xolotlCore::pi;
        if (reg.isSimplex()) {
            Composition comp(reg.getOrigin());
            // Compute the reaction radius
            // TODO: change the hard coded value to get the density from the network/options
            radius = pow(
                (3.0 * (double) comp[Species::Xe]) / (FourPi * _density),
                (1.0 / 3.0));
            if (comp[Species::Xe] == 1)
                radius = impurityRadius;
        }
        else {
            // Loop on the range
            for (auto j : makeIntervalRange(reg[Species::Xe])) {
                radius += pow(
                            (3.0 * (double) j) / (FourPi * _density),
                            (1.0 / 3.0));
            }
            // Average the radius
            radius /= reg.volume();
        }

        return radius;
    }

private:
    AmountType _maxXe;
    double _xeDiffusivity;
    bool _xeDiffusive;
    AmountType _groupingMin;
    AmountType _groupingWidth;
    double _density;
};
}
}
