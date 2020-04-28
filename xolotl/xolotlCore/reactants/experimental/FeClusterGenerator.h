#pragma once

namespace xolotlCore
{
namespace experimental
{

class FeClusterGenerator :
    public
    plsm::refine::Detector<FeClusterGenerator>
{
public:
    using Species = FeSpeciesList;
    using Superclass = plsm::refine::Detector<FeClusterGenerator>;
    using NetworkType = FeReactionNetwork;

    template <typename PlsmContext>
    using Cluster = typename NetworkType::Cluster<PlsmContext>;

    using Region = typename NetworkType::Region;
    using Composition = typename NetworkType::Composition;
    using AmountType = typename NetworkType::AmountType;

    FeClusterGenerator(const IOptions& options)
        :
        _maxHe(options.getMaxImpurity()),
        _maxV(options.getMaxV()),
        _groupingMin(options.getGroupingMin()),
        _groupingWidthHe(options.getGroupingWidthA()),
        _groupingWidthV(options.getGroupingWidthB())
    {
    }

    FeClusterGenerator(const IOptions& options, std::size_t refineDepth)
        :
        Superclass(refineDepth),
        _maxHe(options.getMaxImpurity()),
        _maxV(options.getMaxV()),
        _groupingMin(options.getGroupingMin()),
        _groupingWidthHe(options.getGroupingWidthA()),
        _groupingWidthV(options.getGroupingWidthB())
    {
    }

    KOKKOS_INLINE_FUNCTION
    bool
    intersect(const Region& region) const;

    KOKKOS_INLINE_FUNCTION
    bool
    select(const Region& region) const;

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getFormationEnergy(const Cluster<PlsmContext>& cluster) const noexcept {
        // Always return 0.0 here because we use capillarity laws for the binding energies
        return 0.0;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getMigrationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getDiffusionFactor(const Cluster<PlsmContext>& cluster) const noexcept;

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getReactionRadius(const Cluster<PlsmContext>& cluster,
        double latticeParameter, double interstitialBias, double impurityRadius)
        const noexcept;

private:
    KOKKOS_INLINE_FUNCTION
    double
    getHeVFormationEnergy(Composition comp) const noexcept;

private:
    // Maximum size of single species
    AmountType _maxHe {8};
    AmountType _maxV {0};
    AmountType _groupingMin;
    AmountType _groupingWidthHe;
    AmountType _groupingWidthV;
};
}
}

#include <experimental/FeClusterGenerator.inl>
