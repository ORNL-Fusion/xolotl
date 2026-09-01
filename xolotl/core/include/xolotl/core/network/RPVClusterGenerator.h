#pragma once

#include <plsm/refine/Detector.h>

#include <xolotl/core/network/RPVTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class RPVClusterGenerator : public plsm::refine::Detector<RPVClusterGenerator>
{
public:
	using Species = RPVSpeciesList;
	using Superclass = plsm::refine::Detector<RPVClusterGenerator>;
	using NetworkType = RPVReactionNetwork;

	template <typename PlsmContext>
	using Cluster = typename NetworkType::Cluster<PlsmContext>;

	using Region = typename NetworkType::Region;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;
	using BoolArray = plsm::refine::BoolVec<Region>;

	RPVClusterGenerator(const options::IOptions& options);

	RPVClusterGenerator(
		const options::IOptions& options, std::size_t refineDepth);

	KOKKOS_INLINE_FUNCTION
	bool
	refine(const Region& region, BoolArray& result) const;

	KOKKOS_INLINE_FUNCTION
	bool
	select(const Region& region) const;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getFormationEnergy(const Cluster<PlsmContext>& cluster,
		double latticeParameter, double interstitialBias,
		double impurityRadius) const noexcept
	{
		// Formation energies are defined for solute/precipitates only
		auto clReg = cluster.getRegion();
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();

		constexpr double sEnergy[15] = {0.0, 0.21, 0.37, 0.55, 0.73, 0.92, 1.10,
			1.29, 1.45, 1.60, 1.77, 1.94, 2.12, 2.3, 2.45};

		if (lo.isOnAxis(Species::S)) {
			// Smaller sizes
			if (lo[Species::S] < 15)
				return sEnergy[lo[Species::S]];
			// J to eV
			double conversionFactor = 6.241509e18;
			double sigma = 0.27 * conversionFactor * 1.0e-18; // eV / nm2
			double radius = this->getReactionRadius(
				cluster, latticeParameter, interstitialBias, impurityRadius);
			return 4.0 * ::xolotl::core::pi * radius * radius * sigma;
		}

		return 0.0;
	}

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getMigrationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getDiffusionFactor(const Cluster<PlsmContext>& cluster,
		double latticeParameter) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getReactionRadius(const Cluster<PlsmContext>& cluster,
		double latticeParameter, double interstitialBias,
		double impurityRadius) const noexcept;

private:
	KOKKOS_INLINE_FUNCTION
	double
	getHeVFormationEnergy(Composition comp) const noexcept;

private:
	// Maximum size of single species
	AmountType _maxSize{0};
	AmountType _maxI{0};
	AmountType _groupingMin;
	AmountType _groupingWidth;
};
} // namespace network
} // namespace core
} // namespace xolotl
