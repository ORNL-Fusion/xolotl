#pragma once

#include <plsm/refine/Detector.h>

namespace xolotlCore
{
namespace experimental
{
class AlloyClusterGenerator : public plsm::refine::Detector<AlloyClusterGenerator>
{
public:
    using Superclass = plsm::refine::Detector<AlloyClusterGenerator>;
    using NetworkType = AlloyReactionNetwork;
    using Species = AlloySpecies;

    template <typename PlsmContext>
    using Cluster = Cluster<NetworkType, PlsmContext>;

    using Region = typename NetworkType::Region;
    using Composition = typename NetworkType::Composition;
    using AmountType = typename NetworkType::AmountType;

    AlloyClusterGenerator(const IOptions& options)
        :
        _maxV(options.getMaxV()),
        _maxI(options.getMaxI()),
        _maxSize(options.getMaxImpurity()),
        _groupingMin(options.getGroupingMin()),
        _groupingWidth(options.getGroupingWidthA())
    {
    }

    AlloyClusterGenerator(const IOptions& options, std::size_t refineDepth)
        :
        Superclass(refineDepth),
        _maxV(options.getMaxV()),
        _maxI(options.getMaxI()),
        _maxSize(options.getMaxImpurity()),
        _groupingMin(options.getGroupingMin()),
        _groupingWidth(options.getGroupingWidthA())
    {
    }

    KOKKOS_INLINE_FUNCTION
    bool
    intersect(const Region& region) const
    {
//        if (region[Species::Xe].begin() < _groupingMin) {
//            return true;
//        }
//        if (region[Species::Xe].end() > _maxXe) {
//            return true;
//        }
//        if (region[Species::Xe].length() < max((double) (_groupingWidth + 1), region[Species::Xe].begin() * 1.0e-2)) {
//            return false;
//        }
        return true;
    }

    KOKKOS_INLINE_FUNCTION
    bool
    select(const Region& region) const
    {
        // Each cluster should be on one axis and one axis only
        int nAxis = (region[Species::V].begin() > 0) + (region[Species::I].begin() > 0)
        		+ (region[Species::Perfect].begin() > 0) + (region[Species::Frank].begin() > 0)
				+ (region[Species::Faulted].begin() > 0) + (region[Species::Void].begin() > 0);
        if (nAxis != 1) {
            return false;
        }

        // I
        if (region[Species::I].begin() > _maxI) return false;

        // V
        if (region[Species::V].begin() > _maxV) return false;

        // Perfect
        if (region[Species::Perfect].begin() > 0 && region[Species::Perfect].begin() < _maxI) return false;
        if (region[Species::Perfect].begin() > 0 && region[Species::Perfect].begin() > 45) return false;

        // Frank
        if (region[Species::Frank].begin() > 0 && region[Species::Frank].begin() <= _maxI) return false;
        if (region[Species::Frank].begin() > 0 && region[Species::Frank].begin() > _maxSize) return false;

        // Faulted
        if (region[Species::Faulted].begin() > 0 && region[Species::Faulted].begin() <= _maxV) return false;
        if (region[Species::Faulted].begin() > 0 && region[Species::Faulted].begin() > _maxSize) return false;

        // Void
        if (region[Species::Void].begin() > 0 && region[Species::Void].begin() <= _maxV) return false;
        if (region[Species::Void].begin() > 0 && region[Species::Void].begin() > _maxSize) return false;

        return true;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getFormationEnergy(const Cluster<PlsmContext>& cluster) const noexcept
    {
        const auto& reg = cluster.getRegion();
        Composition comp(reg.getOrigin());
        if (comp.isOnAxis(Species::Perfect))
    		return 3.4 + 2.0 * (pow((double) comp[Species::Perfect], 2.0 / 3.0) - 1.0);
        if (comp.isOnAxis(Species::Frank))
    		return 3.4 + 2.0 * (pow((double) comp[Species::Frank], 2.0 / 3.0) - 1.0);
    	if (comp.isOnAxis(Species::Faulted))
    		return 1.9 + 2.0 * (pow((double) comp[Species::Faulted], 2.0 / 3.0) - 1.0);
    	if (comp.isOnAxis(Species::Void))
    		return 1.9 + 3.4 * (pow((double) comp[Species::Void], 2.0 / 3.0) - 1.0);
    	if (comp.isOnAxis(Species::V))
    		return 1.9 + 3.4 * (pow((double) comp[Species::V], 2.0 / 3.0) - 1.0);
    	if (comp.isOnAxis(Species::I))
    		return 3.4 + 3.5 * (pow((double) comp[Species::I], 2.0 / 3.0) - 1.0);
    	return 0.0;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getMigrationEnergy(const Cluster<PlsmContext>& cluster) const noexcept
    {
        const auto& reg = cluster.getRegion();
        Composition comp(reg.getOrigin());
        double migrationEnergy = std::numeric_limits<double>::infinity();
        if (comp.isOnAxis(Species::Perfect)) {
            return 0.7;
        }
        if (comp.isOnAxis(Species::V)) {
            return 1.3;
        }
        if (comp.isOnAxis(Species::I)) {
            return 0.5;
        }
        return migrationEnergy;
    }

    template <typename PlsmContext>
    KOKKOS_INLINE_FUNCTION
    double
    getDiffusionFactor(const Cluster<PlsmContext>& cluster,
            double latticeParameter) const noexcept
    {
        const auto& reg = cluster.getRegion();
        Composition comp(reg.getOrigin());
        double diffusionFactor = 0.0;
        if (comp.isOnAxis(Species::Perfect)) {
            if (comp[Species::Perfect] < 70) {
                const double jumpDistance = latticeParameter / sqrt(2.0);
                constexpr double phononFrequency = 9.6e12;
                constexpr double jumpsPerPhonon = 1.0;
                constexpr double prefactorExponent = -1.0;
                return phononFrequency * jumpsPerPhonon * jumpDistance
                    * jumpDistance * pow((double) comp[Species::Perfect], prefactorExponent)
                    / (6.0);
            }
            else return diffusionFactor;
        }
        if (comp.isOnAxis(Species::V)) {
            const double jumpDistance = latticeParameter / sqrt(2.0);
            constexpr double phononFrequency = 9.6e12;
            constexpr double jumpsPerPhonon = 1.0;
            constexpr double prefactorExponent = -1.0;
            return phononFrequency * jumpsPerPhonon * jumpDistance
                * jumpDistance * pow((double) comp[Species::V], prefactorExponent)
                / (6.0);
        }
        if (comp.isOnAxis(Species::I)) {
            const double jumpDistance = latticeParameter / sqrt(2.0);
            constexpr double phononFrequency = 9.6e12;
            constexpr double jumpsPerPhonon = 1.0;
            constexpr double prefactorExponent = -1.0;
            return phononFrequency * jumpsPerPhonon * jumpDistance
                * jumpDistance * pow((double) comp[Species::I], prefactorExponent)
                / (6.0);
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
        double radius = 0.0;
        const double prefactor = 0.25 * latticeParameter * latticeParameter
		/ xolotlCore::pi;
        const auto& reg = cluster.getRegion();
        Composition comp(reg.getOrigin());
        if (comp.isOnAxis(Species::Perfect))
    		return sqrt(((double) comp[Species::Perfect] * prefactor) / xolotlCore::perfectBurgers); 0.5 * latticeParameter;
        if (comp.isOnAxis(Species::Frank))
    		return sqrt(((double) comp[Species::Frank] * prefactor) / xolotlCore::frankBurgers); 0.5 * latticeParameter;
    	if (comp.isOnAxis(Species::Faulted))
    		return sqrt(((double) comp[Species::Faulted] * prefactor) / xolotlCore::faultedBurgers); 0.5 * latticeParameter;
    	if (comp.isOnAxis(Species::Void))
    		return pow(0.75 * prefactor * latticeParameter * (double) comp[Species::Void],
    				1.0 / 3.0);
    	if (comp.isOnAxis(Species::V))
    		return pow(0.75 * prefactor * latticeParameter * (double) comp[Species::V],
    				1.0 / 3.0);
    	if (comp.isOnAxis(Species::I))
    		return pow(0.75 * prefactor * latticeParameter * (double) comp[Species::I],
    				1.0 / 3.0);;

        return radius;
    }

private:
    AmountType _maxV;
    AmountType _maxI;
    AmountType _maxSize;
    AmountType _groupingMin;
    AmountType _groupingWidth;
};
}
}
