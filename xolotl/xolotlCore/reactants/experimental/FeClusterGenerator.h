#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class FeClusterGenerator
{
};

template <>
class FeClusterGenerator<FeFullSpeciesList> :
    public
    plsm::refine::Detector<FeClusterGenerator<FeFullSpeciesList>>
{
public:
    using Species = FeFullSpeciesList;
    using Superclass = plsm::refine::Detector<FeClusterGenerator<Species>>;
    using NetworkType = FeReactionNetwork<Species>;

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
        _groupingWidthA(options.getGroupingWidthA()),
        _groupingWidthB(options.getGroupingWidthB())
    {
    }

    FeClusterGenerator(const IOptions& options, std::size_t refineDepth)
        :
        Superclass(refineDepth),
        _maxHe(options.getMaxImpurity()),
        _maxV(options.getMaxV()),
        _groupingMin(options.getGroupingMin()),
        _groupingWidthA(options.getGroupingWidthA()),
        _groupingWidthB(options.getGroupingWidthB())
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
    AmountType _groupingWidthA;
    AmountType _groupingWidthB;
};
}
}

#include <experimental/FeClusterGenerator.inl>
