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
IReactionNetwork::AmountType
getMaxHePerV(IReactionNetwork::AmountType amtV, double ratio) noexcept
{
	using AmountType = IReactionNetwork::AmountType;

	/**
	 * The maximum number of helium atoms that can be combined with a
	 * vacancy cluster with size equal to the index i.
	 * It could support a mixture of up to nine
	 * helium atoms with one vacancy.
	 */
	constexpr Kokkos::Array<AmountType, 30> maxHePerV = {0, 9, 14, 18, 20, 27,
		30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 98, 100, 101,
		103, 105, 107, 109, 110, 112, 116};

	if (amtV < maxHePerV.size()) {
		return maxHePerV[amtV];
	}
	return util::max((AmountType)(ratio * amtV),
		maxHePerV[maxHePerV.size() - 1] + amtV - (AmountType)maxHePerV.size() +
			1);
}

KOKKOS_INLINE_FUNCTION
double
getMaxHPerV(double amtV, double latticeParameter, double temp) noexcept
{
	// Special case for 1 V
	if (amtV < 1.5)
		return 6;

	// Compute the radius first (in nm)
	double rB = (sqrt(3.0) / 4.0) * latticeParameter +
		pow((3.0 * pow(latticeParameter, 3.0) * (double)amtV) /
				(8.0 * ::xolotl::core::pi),
			(1.0 / 3.0)) -
		pow((3.0 * pow(latticeParameter, 3.0)) / (8.0 * ::xolotl::core::pi),
			(1.0 / 3.0));

	// Radius in m
	double rBm = rB * 1.0e-9;

	// Compute the shear modulus (in Pa)
	double G = 163.4e9 * (1.0 - 0.18 * (temp / 3700.0));

	// Compute the loop punching pressure (in MPa)
	double p = ((2.0 * 2.65 / rBm) + (G * 3.0e-10 / rBm)) * 1.0e-6;

	// Equation of state to get the molar volume (in cm3 / mol)
	double v = 176.33 * pow(p, -1.0 / 3.0) - 633.675 * pow(p, -2.0 / 3.0) -
		304.574 * pow(p, -4.0 / 3.0) + (731.393 + 8.59805 * temp) / p;

	// Get the number of hydrogen molecules
	double nM = 6.02214e23 * (4.0 * ::xolotl::core::pi * rBm * rBm * rBm) /
		(3.0 * v * 1.0e-6);

	return 2.0 * nM;
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

	// KOKKOS_FUNCTION
	// static AmountType
	// getMaxHePerV(AmountType amtV, double ratio) noexcept;

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
	AmountType _maxPureV{0};
	AmountType _maxI{0};
	AmountType _groupingMin;
	AmountType _groupingWidthA;
	AmountType _groupingWidthB;
	double _hevRatio{4.0};

	// The temperature
	double _temperature{933.0};
	// The lattice parameter
	double _lattice{xolotl::core::tungstenLatticeConstant};
};
} // namespace network
} // namespace core
} // namespace xolotl
