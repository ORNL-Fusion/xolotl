#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIClusterGenerator
{
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

    template <typename PlsmContext>
    using Cluster = typename NetworkType::Cluster<PlsmContext>;

    using Region = typename NetworkType::Region;
    using Composition = typename NetworkType::Composition;
    using AmountType = typename NetworkType::AmountType;

    PSIClusterGenerator(const IOptions& options)
        :
        _hydrogenRadiusFactor(options.getHydrogenFactor()),
        _maxHe(options.getMaxImpurity()),
        _maxD(options.getMaxD()),
        _maxT(options.getMaxT()),
        _maxV(options.getMaxV()),
        _groupingMin(options.getGroupingMin()),
        _groupingWidthA(options.getGroupingWidthA()),
        _groupingWidthB(options.getGroupingWidthB())
    {
    }

    PSIClusterGenerator(const IOptions& options, std::size_t refineDepth)
        :
        Superclass(refineDepth),
        _hydrogenRadiusFactor(options.getHydrogenFactor()),
        _maxHe(options.getMaxImpurity()),
        _maxD(options.getMaxD()),
        _maxT(options.getMaxT()),
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

    KOKKOS_INLINE_FUNCTION
    static
    AmountType
    getMaxHePerV(AmountType amtV) noexcept;

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getFormationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

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
    // The factor between He and H radius sizes
    double _hydrogenRadiusFactor {0.25};

    // Maximum size of single species
    AmountType _maxHe {8};
    AmountType _maxD {1};
    AmountType _maxT {1};
    AmountType _maxV {0};
    AmountType _groupingMin;
    AmountType _groupingWidthA;
    AmountType _groupingWidthB;
};
}
}

#include <experimental/PSIClusterGenerator.inl>
