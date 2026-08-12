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

/**
 * Continuous maximum number of He per V, for use in the single size bubble
 * model (SSBM).
 *
 * IMPORTANT: this is deliberately distinct from getMaxHePerV() above. That
 * function takes an INTEGER vacancy count and returns an INTEGER from a lookup
 * table, which is correct for generating the discrete cluster network. The SSBM
 * evaluates the He/V closure at the *average* V of the large bubble, which is a
 * continuous solution variable (moment / concentration). Truncating it to an
 * integer makes the residual a staircase in the state vector: the sigmoid
 * centre jumps every time <V> crosses an integer, the Jacobian is zero almost
 * everywhere and undefined on the jumps, and Newton has nothing to converge to.
 *
 * This function must therefore stay C^1 in amtV, and getMaxHePerVdV() below
 * must remain its exact derivative.
 *
 * TODO: the current form reproduces the asymptotic ratio only. Replace the
 * body with a proper He equation of state -- Hammond et al., Sci. Rep. 10
 * (2020) 2192 is parameter-free and was written for coarse-grained models like
 * this one; the MLB-EOS on feature-PSI-large-bubble is the other candidate.
 * When you do, derive getMaxHePerVdV analytically the same way getMaxHPerVdV
 * is derived (chain rule through r_B -> p -> v_molar -> n_molecules).
 *
 * @param amtV The average number of vacancies in the bubble (continuous)
 * @param latticeParameter The lattice parameter
 * @param temp The temperature
 * @return The maximum number of He the bubble can hold at this size
 */
KOKKOS_INLINE_FUNCTION
double
getMaxHePerVCont(double amtV, double latticeParameter, double temp) noexcept
{
	// Floor. Newton can produce slightly negative <V> on intermediate
	// iterates; clamping here keeps the function defined and keeps the
	// derivative consistent with getMaxHePerVdV.
	//
	// KNOWN LIMITATION: this clamp makes the function C^0 but not C^1 -- there
	// is a slope discontinuity exactly at amtV == 1.0 (below it the value is
	// constant, above it linear). Verified against central differences: the
	// analytic derivative agrees to machine precision everywhere except at
	// amtV == 1.0, where the one-sided FD is half the analytic value.
	//
	// This is enormously better than truncating to the integer lookup table,
	// which is discontinuous at EVERY integer, and a kink on a measure-zero
	// set is normally harmless for Newton. But if you see the solver stalling
	// with <V> pinned near 1, smooth the floor, e.g.
	//     double v = 0.5 * (amtV + sqrt(amtV * amtV + eps));   // softplus-like
	// and update getMaxHePerVdV to match.
	double v = util::max(amtV, 1.0);

	constexpr double hevRatio = 4.0;
	return hevRatio * v;
}

/**
 * Derivative of getMaxHePerVCont with respect to amtV.
 * Must be kept exactly consistent with getMaxHePerVCont: the Jacobian entries
 * for the SSBM trap mutation sigmoid chain-rule through this.
 */
KOKKOS_INLINE_FUNCTION
double
getMaxHePerVdV(double amtV, double latticeParameter, double temp) noexcept
{
	// Below the floor applied in getMaxHePerVCont the function is constant
	if (amtV < 1.0)
		return 0.0;

	constexpr double hevRatio = 4.0;
	return hevRatio;
}

KOKKOS_INLINE_FUNCTION
double
getMaxHPerV(double amtV, double latticeParameter, double temp) noexcept
{
	// Special case for 1 V
	if (amtV < 1.5 and amtV >= 0.0)
		return 6.0;

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

KOKKOS_INLINE_FUNCTION
double
getMaxHPerVdV(double amtV, double latticeParameter, double temp) noexcept
{
	// Special case for 1 V
	if (amtV < 1.5 and amtV >= 0.0)
		return 0.0;

	// Compute the radius first (in nm)
	double rB = (sqrt(3.0) / 4.0) * latticeParameter +
		pow((3.0 * pow(latticeParameter, 3.0) * (double)amtV) /
				(8.0 * ::xolotl::core::pi),
			(1.0 / 3.0)) -
		pow((3.0 * pow(latticeParameter, 3.0)) / (8.0 * ::xolotl::core::pi),
			(1.0 / 3.0));

	// Radius in m
	double rBm = rB * 1.0e-9;
	double rBmdV = 1.0e-9 *
		(pow(latticeParameter, 3.0) / (8.0 * ::xolotl::core::pi)) *
		pow((3.0 * pow(latticeParameter, 3.0) * (double)amtV) /
				(8.0 * ::xolotl::core::pi),
			(-2.0 / 3.0));

	// Compute the shear modulus (in Pa)
	double G = 163.4e9 * (1.0 - 0.18 * (temp / 3700.0));

	// Compute the loop punching pressure (in MPa)
	double p = ((2.0 * 2.65 / rBm) + (G * 3.0e-10 / rBm)) * 1.0e-6;
	double pdV = -(2.0 * 2.65 + G * 3.0e-10) * 1.0e-6 * rBmdV / (rBm * rBm);

	// Equation of state to get the molar volume (in cm3 / mol)
	double v = 176.33 * pow(p, -1.0 / 3.0) - 633.675 * pow(p, -2.0 / 3.0) -
		304.574 * pow(p, -4.0 / 3.0) + (731.393 + 8.59805 * temp) / p;
	double vdV = -(176.33 / 3.0) * pdV * pow(p, -4.0 / 3.0) +
		((2.0 * 633.675) / 3.0) * pdV * pow(p, -5.0 / 3.0) +
		((4.0 * 304.574) / 3.0) * pdV * pow(p, -7.0 / 3.0) -
		(731.393 + 8.59805 * temp) * pdV / (p * p);

	// Get the number of hydrogen molecules
	double m = rBm * rBm * rBm;
	double mdV = 3.0 * rBmdV * rBm * rBm;
	double nMdV = (6.02214e23 * (4.0 * ::xolotl::core::pi) / (3.0 * 1.0e-6)) *
		(mdV * v - m * vdV) / (v * v);

	return 2.0 * nMdV;
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
