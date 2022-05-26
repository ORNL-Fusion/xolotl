#pragma once

#include <xolotl/core/Constants.h>
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
getMaxHePerV(IReactionNetwork::AmountType amtV) noexcept
{
	using AmountType = IReactionNetwork::AmountType;

	/**
	 * The maximum number of helium atoms that can be combined with a
	 * vacancy cluster with size equal to the index i.
	 * It could support a mixture of up to nine
	 * helium atoms with one vacancy.
	 */
	constexpr Kokkos::Array<AmountType, 15> maxHePerV = {
		0, 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70};

	if (amtV < maxHePerV.size()) {
		return maxHePerV[amtV];
	}

	return util::max((AmountType)(pow((double)amtV, 0.86) * 5.0),
		maxHePerV[maxHePerV.size() - 1] + amtV - (AmountType)maxHePerV.size() +
			1);
}

KOKKOS_INLINE_FUNCTION
double
computeBenedictF1(const double temp) noexcept
{
	return xolotl::core::A_11 / sqrt(temp) + xolotl::core::A01 +
		xolotl::core::A21 * temp;
}

KOKKOS_INLINE_FUNCTION
double
computeBenedictF2(const double temp) noexcept
{
	return xolotl::core::A02 + xolotl::core::A22 * temp;
}

KOKKOS_INLINE_FUNCTION
double
computeBenedictF3(const double temp) noexcept
{
	return xolotl::core::A_23 / temp + xolotl::core::A_13 / sqrt(temp) +
		xolotl::core::A03 + xolotl::core::A23 * temp;
}

KOKKOS_INLINE_FUNCTION
IReactionNetwork::AmountType
getMaxHePerVEq(IReactionNetwork::AmountType amtV, double latticeParameter,
	double temp) noexcept
{
	using AmountType = IReactionNetwork::AmountType;
	constexpr Kokkos::Array<AmountType, 15> maxHePerV = {
		0, 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70};

	if (amtV < maxHePerV.size()) {
		return maxHePerV[amtV];
	}

	double omega =
		0.5 * latticeParameter * latticeParameter * latticeParameter * 1.0e-27;
	double nume = amtV * omega;
	double Fs = pow(3.0 * omega / (4.0 * xolotl::core::pi), 1.0 / 3.0);
	double term1 = computeBenedictF1(temp) * pow(Fs, 1.0 / 3.0) *
		pow((double)amtV, 1.0 / 9.0) *
		pow(2.0 * xolotl::core::gammaEOS, -1.0 / 3.0);
	double term2 = computeBenedictF2(temp) * pow(Fs, 2.0 / 3.0) *
		pow((double)amtV, 2.0 / 9.0) *
		pow(2.0 * xolotl::core::gammaEOS, -2.0 / 3.0);
	double term3 = computeBenedictF3(temp) * Fs * pow((double)amtV, 1.0 / 3.0) /
		(2.0 * xolotl::core::gammaEOS);

	return util::max((AmountType)(nume / (term1 + term2 + term3)),
		maxHePerV[maxHePerV.size() - 1] + amtV - (AmountType)maxHePerV.size() +
			1);
}

KOKKOS_INLINE_FUNCTION
IReactionNetwork::AmountType
getMaxHePerVLoop(IReactionNetwork::AmountType amtV, double latticeParameter,
	double temp) noexcept
{
	using AmountType = IReactionNetwork::AmountType;
	constexpr Kokkos::Array<AmountType, 15> maxHePerV = {
		0, 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70};

	if (amtV < maxHePerV.size()) {
		return maxHePerV[amtV];
	}

	double omega =
		0.5 * latticeParameter * latticeParameter * latticeParameter * 1.0e-27;
	double bEOS = latticeParameter * 0.5 * sqrt(3.0) * 1.0e-9;
	double nume = amtV * omega;
	double Fs = pow(3.0 * omega / (4.0 * xolotl::core::pi), 1.0 / 3.0);
	double term1 = computeBenedictF1(temp) * pow(Fs, 1.0 / 3.0) *
		pow((double)amtV, 1.0 / 9.0) *
		pow(2.0 * xolotl::core::gammaEOS + xolotl::core::gEOS * bEOS,
			-1.0 / 3.0);
	double term2 = computeBenedictF2(temp) * pow(Fs, 2.0 / 3.0) *
		pow((double)amtV, 2.0 / 9.0) *
		pow(2.0 * xolotl::core::gammaEOS + xolotl::core::gEOS * bEOS,
			-2.0 / 3.0);
	double term3 = computeBenedictF3(temp) * Fs * pow((double)amtV, 1.0 / 3.0) /
		(2.0 * xolotl::core::gammaEOS + xolotl::core::gEOS * bEOS);

	return util::max((AmountType)(nume / (term1 + term2 + term3)),
		maxHePerV[maxHePerV.size() - 1] + amtV - (AmountType)maxHePerV.size() +
			1);
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
	using BoolArray = typename Superclass::template BoolVec<Region>;

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

	// The temperature
	double _temperature{933.0};
	// The lattice parameter
	double _lattice{xolotl::core::tungstenLatticeConstant};
};
} // namespace network
} // namespace core
} // namespace xolotl
