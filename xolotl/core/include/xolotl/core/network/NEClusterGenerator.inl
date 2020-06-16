#pragma once

#include <xolotl/util/MathUtils.h>
#include <xolotl/core/Constants.h>

namespace xolotl
{
namespace core
{
namespace network
{
NEClusterGenerator::NEClusterGenerator(const options::IOptions& opts)
    :
    _maxXe(opts.getMaxImpurity()),
    _xeDiffusivity(opts.getXenonDiffusivity()),
    _xeDiffusive(_xeDiffusivity > 0.0),
    _groupingMin(opts.getGroupingMin()),
    _groupingWidth(opts.getGroupingWidthA()),
    _density(opts.getDensity())
{
}

NEClusterGenerator::NEClusterGenerator(const options::IOptions& opts,
        std::size_t refineDepth)
    :
    Superclass(refineDepth),
    _maxXe(opts.getMaxImpurity()),
    _xeDiffusivity(opts.getXenonDiffusivity()),
    _xeDiffusive(_xeDiffusivity > 0.0),
    _groupingMin(opts.getGroupingMin()),
    _groupingWidth(opts.getGroupingWidthA()),
    _density(opts.getDensity())
{
}

KOKKOS_INLINE_FUNCTION
bool
NEClusterGenerator::intersect(const Region& region) const
{
    if (region[Species::Xe].begin() < _groupingMin) {
        return true;
    }
    if (region[Species::Xe].end() > _maxXe) {
        return true;
    }
    if (region[Species::Xe].length() < util::max((double) (_groupingWidth + 1), region[Species::Xe].begin() * 1.0e-2)) {
        return false;
    }
    return true;
}

KOKKOS_INLINE_FUNCTION
bool
NEClusterGenerator::select(const Region& region) const
{
    // Remove 0
    if (region[Species::Xe].end() == 1) {
        return false;
    }

    if (region[Species::Xe].begin() > _maxXe) {
        return false;
    }

    return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
NEClusterGenerator::getFormationEnergy(const Cluster<PlsmContext>& cluster)
    const noexcept
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
    return 79.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
NEClusterGenerator::getMigrationEnergy(const Cluster<PlsmContext>& cluster)
    const noexcept
{
    constexpr double xeOneMigration = 1.0;

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
NEClusterGenerator::getDiffusionFactor(const Cluster<PlsmContext>& cluster,
    double latticeParameter) const noexcept
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
NEClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
    double latticeParameter, double interstitialBias, double impurityRadius)
    const noexcept
{
    const auto& reg = cluster.getRegion();
    double radius = 0.0;
    double FourPi = 4.0 * ::xolotl::core::pi;
    if (reg.isSimplex()) {
        Composition comp(reg.getOrigin());
        // Compute the reaction radius
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
}
}
}
