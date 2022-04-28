#pragma once

#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace psi
{
KOKKOS_INLINE_FUNCTION
double
getMaxHePerV(double amtV) noexcept
{
	using AmountType = IReactionNetwork::AmountType;

	amtV = std::max(0.0, amtV);

	/**
	 * The maximum number of helium atoms that can be combined with a
	 * vacancy cluster with size equal to the index i.
	 * It could support a mixture of up to nine
	 * helium atoms with one vacancy.
	 */
	constexpr Kokkos::Array<AmountType, 15> maxHePerV = {
		8, 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70};

	if (amtV < maxHePerV.size()) {
		return maxHePerV[(AmountType)amtV];
	}

	return util::max(pow(amtV, 0.86) * 5.0,
		maxHePerV[maxHePerV.size() - 1] + amtV - maxHePerV.size() + 1);
}
} // namespace psi

template <typename TSpeciesEnum>
class PSIClusterGenerator :
	public plsm::refine::Detector<PSIClusterGenerator<TSpeciesEnum>>
{
public:
	using Species = TSpeciesEnum;
	using Superclass = plsm::refine::Detector<PSIClusterGenerator<Species>>;
	using NetworkType = PSIReactionNetwork<Species>;

	template <typename PlsmContext>
	using Cluster = typename NetworkType::template Cluster<PlsmContext>;

	using Region = typename NetworkType::Region;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;
	using BoolArray = plsm::refine::BoolVec<Region>;

	PSIClusterGenerator(const options::IOptions& opts);

	PSIClusterGenerator(const options::IOptions& opts, std::size_t refineDepth);

	KOKKOS_INLINE_FUNCTION
	bool
	refine(const Region& region, BoolArray& result) const;

	KOKKOS_INLINE_FUNCTION
	bool
	select(const Region& region) const;

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
	getDiffusionFactor(const Cluster<PlsmContext>& cluster,
		double latticeParameter) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getReactionRadius(const Cluster<PlsmContext>& cluster,
		double latticeParameter, double interstitialBias,
		double impurityRadius) const noexcept;

	KOKKOS_INLINE_FUNCTION
	static double
	getHeVFormationEnergy(Composition comp);

private:
	// The factor between He and H radius sizes
	double _hydrogenRadiusFactor{0.25};

	// Maximum size of single species
	AmountType _maxHe{8};
	AmountType _maxD{1};
	AmountType _maxT{1};
	AmountType _maxV{0};
	AmountType _maxI{0};
	AmountType _groupingMin;
	AmountType _groupingWidthA;
	AmountType _groupingWidthB;
};
} // namespace network
} // namespace core
} // namespace xolotl
