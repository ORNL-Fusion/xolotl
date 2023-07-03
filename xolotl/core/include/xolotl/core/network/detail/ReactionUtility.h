#pragma once

#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief Reflected region are used to combine vacancy and interstitial spaces
 * as they are opposite.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
Kokkos::Array<
	plsm::Region<plsm::DifferenceType<typename TRegion::ScalarType>, Dim>, 4>
initReflectedRegions(const TRegion& cl1Reg, const TRegion& cl2Reg,
	const TRegion& pr1Reg, const TRegion& pr2Reg)
{
	using Species = typename TRegion::EnumIndex;
	using SpeciesSequence = SpeciesEnumSequence<Species, TRegion::dimension()>;
	using SpeciesRange = EnumSequenceRange<Species, TRegion::dimension()>;
	using Ival =
		plsm::Interval<plsm::DifferenceType<typename TRegion::ScalarType>>;

	constexpr auto speciesRange =
		SpeciesRange(SpeciesSequence::first(), SpeciesSequence::lastNoI());
	plsm::Region<plsm::DifferenceType<typename TRegion::ScalarType>, Dim> cl1RR,
		cl2RR, pr1RR, pr2RR;
	for (auto i : speciesRange) {
		cl1RR[i()] = Ival(cl1Reg[i].begin(), cl1Reg[i].end());
		cl2RR[i()] = Ival(cl2Reg[i].begin(), cl2Reg[i].end());
		pr1RR[i()] = Ival(pr1Reg[i].begin(), pr1Reg[i].end());
		pr2RR[i()] = Ival(pr2Reg[i].begin(), pr2Reg[i].end());
	}
	return {cl1RR, cl2RR, pr1RR, pr2RR};
}

/**
 * @brief Specific case with more than one type of V. Here it can be reflected
 * of V or I depending on the product type.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() > 2),
	Kokkos::Array<
		plsm::Region<plsm::DifferenceType<typename TRegion::ScalarType>, Dim>,
		4>>
updateReflectedRegionsForCoefs(const TRegion& cl1Reg, const TRegion& cl2Reg,
	const TRegion& pr1Reg, const TRegion& pr2Reg,
	plsm::DifferenceType<typename TRegion::ScalarType> m = 1)
{
	using Species = typename TRegion::EnumIndex;
	using Ival =
		plsm::Interval<plsm::DifferenceType<typename TRegion::ScalarType>>;

	// Initialize the reflected regions
	auto rRegions = initReflectedRegions<Dim>(cl1Reg, cl2Reg, pr1Reg, pr2Reg);
	auto cl1RR = rRegions[0];
	auto cl2RR = rRegions[1];
	auto pr1RR = rRegions[2];
	auto pr2RR = rRegions[3];
	auto vIndex = static_cast<std::underlying_type_t<Species>>(Species::V);
	// The first product tells us how to project
	if (pr1Reg[Species::I].end() > 1 || pr1Reg[Species::Perfect].end() > 1 ||
		pr1Reg[Species::Frank].end() > 1) {
		// Project on I
		cl1RR[vIndex] =
			Ival(cl1Reg[Species::I].begin() + cl1Reg[Species::Perfect].begin() +
					cl1Reg[Species::Frank].begin() - cl1Reg[Species::V].end() -
					cl1Reg[Species::Void].end() -
					cl1Reg[Species::Faulted].end() + 3,
				cl1Reg[Species::I].end() + cl1Reg[Species::Perfect].end() +
					cl1Reg[Species::Frank].end() - cl1Reg[Species::V].begin() -
					cl1Reg[Species::Void].begin() -
					cl1Reg[Species::Faulted].begin() - 2);
		cl2RR[vIndex] =
			Ival(cl2Reg[Species::I].begin() + cl2Reg[Species::Perfect].begin() +
					cl2Reg[Species::Frank].begin() - cl2Reg[Species::V].end() -
					cl2Reg[Species::Void].end() -
					cl2Reg[Species::Faulted].end() + 3,
				cl2Reg[Species::I].end() + cl2Reg[Species::Perfect].end() +
					cl2Reg[Species::Frank].end() - cl2Reg[Species::V].begin() -
					cl2Reg[Species::Void].begin() -
					cl2Reg[Species::Faulted].begin() - 2);
		pr1RR[vIndex] =
			Ival(pr1Reg[Species::I].begin() + pr1Reg[Species::Perfect].begin() +
					pr1Reg[Species::Frank].begin() - pr1Reg[Species::V].end() -
					pr1Reg[Species::Void].end() -
					pr1Reg[Species::Faulted].end() + 3,
				pr1Reg[Species::I].end() + pr1Reg[Species::Perfect].end() +
					pr1Reg[Species::Frank].end() - pr1Reg[Species::V].begin() -
					pr1Reg[Species::Void].begin() -
					pr1Reg[Species::Faulted].begin() - 2);
		pr2RR[vIndex] =
			Ival(pr2Reg[Species::I].begin() + pr2Reg[Species::Perfect].begin() +
					pr2Reg[Species::Frank].begin() - pr2Reg[Species::V].end() -
					pr2Reg[Species::Void].end() -
					pr2Reg[Species::Faulted].end() + 3,
				pr2Reg[Species::I].end() + pr2Reg[Species::Perfect].end() +
					pr2Reg[Species::Frank].end() - pr2Reg[Species::V].begin() -
					pr2Reg[Species::Void].begin() -
					pr2Reg[Species::Faulted].begin() - 2);
	}
	else {
		// Project on V
		cl1RR[vIndex] = Ival(cl1Reg[Species::V].begin() +
				cl1Reg[Species::Void].begin() +
				cl1Reg[Species::Faulted].begin() - cl1Reg[Species::I].end() -
				cl1Reg[Species::Perfect].end() - cl1Reg[Species::Frank].end() +
				3,
			cl1Reg[Species::V].end() + cl1Reg[Species::Void].end() +
				cl1Reg[Species::Faulted].end() - cl1Reg[Species::I].begin() -
				cl1Reg[Species::Perfect].begin() -
				cl1Reg[Species::Frank].begin() - 2);
		cl2RR[vIndex] = Ival(cl2Reg[Species::V].begin() +
				cl2Reg[Species::Void].begin() +
				cl2Reg[Species::Faulted].begin() - cl2Reg[Species::I].end() -
				cl2Reg[Species::Perfect].end() - cl2Reg[Species::Frank].end() +
				3,
			cl2Reg[Species::V].end() + cl2Reg[Species::Void].end() +
				cl2Reg[Species::Faulted].end() - cl2Reg[Species::I].begin() -
				cl2Reg[Species::Perfect].begin() -
				cl2Reg[Species::Frank].begin() - 2);
		pr1RR[vIndex] = Ival(pr1Reg[Species::V].begin() +
				pr1Reg[Species::Void].begin() +
				pr1Reg[Species::Faulted].begin() - pr1Reg[Species::I].end() -
				pr1Reg[Species::Perfect].end() - pr1Reg[Species::Frank].end() +
				3,
			pr1Reg[Species::V].end() + pr1Reg[Species::Void].end() +
				pr1Reg[Species::Faulted].end() - pr1Reg[Species::I].begin() -
				pr1Reg[Species::Perfect].begin() -
				pr1Reg[Species::Frank].begin() - 2);
		pr2RR[vIndex] = Ival(pr2Reg[Species::V].begin() +
				pr2Reg[Species::Void].begin() +
				pr2Reg[Species::Faulted].begin() - pr2Reg[Species::I].end() -
				pr2Reg[Species::Perfect].end() - pr2Reg[Species::Frank].end() +
				3,
			pr2Reg[Species::V].end() + pr2Reg[Species::Void].end() +
				pr2Reg[Species::Faulted].end() - pr2Reg[Species::I].begin() -
				pr2Reg[Species::Perfect].begin() -
				pr2Reg[Species::Frank].begin() - 2);
	}
	return {cl1RR, cl2RR, pr1RR, pr2RR};
}

// JMR edit:

/**
 * @brief Specific case with two typse of V  (Zr - basal). Here it can be
 * reflected of V or I depending on the product type.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() == 2),
	Kokkos::Array<
		plsm::Region<plsm::DifferenceType<typename TRegion::ScalarType>, Dim>,
		4>>
updateReflectedRegionsForCoefs(const TRegion& cl1Reg, const TRegion& cl2Reg,
	const TRegion& pr1Reg, const TRegion& pr2Reg,
	plsm::DifferenceType<typename TRegion::ScalarType> m = 1)
{
	using Species = typename TRegion::EnumIndex;
	using Ival =
		plsm::Interval<plsm::DifferenceType<typename TRegion::ScalarType>>;

	// Initialize the reflected regions
	auto rRegions = initReflectedRegions<Dim>(cl1Reg, cl2Reg, pr1Reg, pr2Reg);
	auto cl1RR = rRegions[0];
	auto cl2RR = rRegions[1];
	auto pr1RR = rRegions[2];
	auto pr2RR = rRegions[3];
	auto vIndex = static_cast<std::underlying_type_t<Species>>(Species::V);
	// The first product tells us how to project
	if (pr1Reg[Species::I].end() > 1) {
		// Project on I
		cl1RR[vIndex] = Ival(cl1Reg[Species::I].begin() -
				cl1Reg[Species::V].end() - cl1Reg[Species::Basal].end() + 2,
			cl1Reg[Species::I].end() - cl1Reg[Species::V].begin() -
				cl1Reg[Species::Basal].begin());
		cl2RR[vIndex] = Ival(cl2Reg[Species::I].begin() -
				cl2Reg[Species::V].end() - cl2Reg[Species::Basal].end() + 2,
			cl2Reg[Species::I].end() - cl2Reg[Species::V].begin() -
				cl2Reg[Species::Basal].begin());
		pr1RR[vIndex] = Ival(pr1Reg[Species::I].begin() -
				pr1Reg[Species::V].end() - pr1Reg[Species::Basal].end() + 2,
			pr1Reg[Species::I].end() - pr1Reg[Species::V].begin() -
				pr1Reg[Species::Basal].begin());
		pr2RR[vIndex] = Ival(pr2Reg[Species::I].begin() -
				pr2Reg[Species::V].end() - pr2Reg[Species::Basal].end() + 2,
			pr2Reg[Species::I].end() - pr2Reg[Species::V].begin() -
				pr2Reg[Species::Basal].begin());
	}
	else {
		// Project on V
		cl1RR[vIndex] = Ival(cl1Reg[Species::V].begin() +
				cl1Reg[Species::Basal].begin() - cl1Reg[Species::I].end() + 1,
			cl1Reg[Species::V].end() + cl1Reg[Species::Basal].end() -
				cl1Reg[Species::I].begin() - 1);
		cl2RR[vIndex] = Ival(cl2Reg[Species::V].begin() +
				cl2Reg[Species::Basal].begin() - cl2Reg[Species::I].end() + 1,
			cl2Reg[Species::V].end() + cl2Reg[Species::Basal].end() -
				cl2Reg[Species::I].begin() - 1);
		pr1RR[vIndex] = Ival(pr1Reg[Species::V].begin() +
				pr1Reg[Species::Basal].begin() - pr1Reg[Species::I].end() + 1,
			pr1Reg[Species::V].end() + pr1Reg[Species::Basal].end() -
				pr1Reg[Species::I].begin() - 1);
		pr2RR[vIndex] = Ival(pr2Reg[Species::V].begin() +
				pr2Reg[Species::Basal].begin() - pr2Reg[Species::I].end() + 1,
			pr2Reg[Species::V].end() + pr2Reg[Species::Basal].end() -
				pr2Reg[Species::I].begin() - 1);
	}
	return {cl1RR, cl2RR, pr1RR, pr2RR};
}

// end JMR edit

/**
 * @brief Specific case with one type of V. Here it can be reflected
 * of V or I depending on the product type.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() == 1),
	Kokkos::Array<
		plsm::Region<plsm::DifferenceType<typename TRegion::ScalarType>, Dim>,
		4>>
updateReflectedRegionsForCoefs(const TRegion& cl1Reg, const TRegion& cl2Reg,
	const TRegion& pr1Reg, const TRegion& pr2Reg,
	plsm::DifferenceType<typename TRegion::ScalarType> m = 1)
{
	using Species = typename TRegion::EnumIndex;
	using Ival =
		plsm::Interval<plsm::DifferenceType<typename TRegion::ScalarType>>;

	// Initialize the reflected regions
	auto rRegions = initReflectedRegions<Dim>(cl1Reg, cl2Reg, pr1Reg, pr2Reg);
	auto cl1RR = rRegions[0];
	auto cl2RR = rRegions[1];
	auto pr1RR = rRegions[2];
	auto pr2RR = rRegions[3];
	auto vIndex = static_cast<std::underlying_type_t<Species>>(Species::V);
	// The first product tells us how to project
	if (pr1Reg[Species::I].end() > 1) {
		// Project on I
		cl1RR[vIndex] =
			Ival(cl1Reg[Species::I].begin() - cl1Reg[Species::V].end() + 1,
				cl1Reg[Species::I].end() - cl1Reg[Species::V].begin());
		cl2RR[vIndex] =
			Ival(cl2Reg[Species::I].begin() - cl2Reg[Species::V].end() + 1,
				cl2Reg[Species::I].end() - cl2Reg[Species::V].begin());
		pr1RR[vIndex] =
			Ival(pr1Reg[Species::I].begin() - pr1Reg[Species::V].end() + 1,
				pr1Reg[Species::I].end() - pr1Reg[Species::V].begin());
		pr2RR[vIndex] =
			Ival(pr2Reg[Species::I].begin() * m - pr2Reg[Species::V].end() + 1,
				(pr2Reg[Species::I].end() - 1) * m + 1 -
					pr2Reg[Species::V].begin());
	}
	else {
		// Project on V
		cl1RR[vIndex] =
			Ival(cl1Reg[Species::V].begin() - cl1Reg[Species::I].end() + 1,
				cl1Reg[Species::V].end() - cl1Reg[Species::I].begin());
		cl2RR[vIndex] =
			Ival(cl2Reg[Species::V].begin() - cl2Reg[Species::I].end() + 1,
				cl2Reg[Species::V].end() - cl2Reg[Species::I].begin());
		pr1RR[vIndex] =
			Ival(pr1Reg[Species::V].begin() - pr1Reg[Species::I].end() + 1,
				pr1Reg[Species::V].end() - pr1Reg[Species::I].begin());
		pr2RR[vIndex] = Ival(
			pr2Reg[Species::V].begin() - (pr2Reg[Species::I].end() - 1) * m,
			pr2Reg[Species::V].end() - (pr2Reg[Species::I].begin() * m));
	}
	return {cl1RR, cl2RR, pr1RR, pr2RR};
}

/**
 * @brief Specific case where there is no V. Simply initialize
 * the reflected regions.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() == 0),
	Kokkos::Array<
		plsm::Region<plsm::DifferenceType<typename TRegion::ScalarType>, Dim>,
		4>>
updateReflectedRegionsForCoefs(const TRegion& cl1Reg, const TRegion& cl2Reg,
	const TRegion& pr1Reg, const TRegion& pr2Reg,
	plsm::DifferenceType<typename TRegion::ScalarType> m = 1)
{
	return initReflectedRegions<Dim>(cl1Reg, cl2Reg, pr1Reg, pr2Reg);
}

/**
 * @brief Specific case with more than one type of V. The I dispersion has to be
 * taken into account.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<((numberOfVacancySpecies<typename TRegion::EnumIndex>() > 1)),
	plsm::SpaceVector<double, Dim>>
getReflectedDispersionForCoefs(const TRegion& clReg)
{
	using Species = typename TRegion::EnumIndex;
	auto disp = clReg.dispersion();
	auto vIndex = static_cast<std::underlying_type_t<Species>>(Species::V);
	for (auto i = 0; i < Dim; i++) {
		if (i == vIndex)
			continue;
		disp[vIndex] += disp[i];
	}
	disp[vIndex] -= Dim - 1.0;
	return disp;
}

/**
 * @brief Specific case with one V. I dispersion is added to V.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() == 1),
	plsm::SpaceVector<double, Dim>>
getReflectedDispersionForCoefs(const TRegion& clReg)
{
	using Species = typename TRegion::EnumIndex;
	auto disp = clReg.dispersion();
	auto vIndex = static_cast<std::underlying_type_t<Species>>(Species::V);
	auto iIndex = static_cast<std::underlying_type_t<Species>>(Species::I);
	disp[vIndex] += disp[iIndex] - 1.0;
	return disp;
}

/**
 * @brief Specific case no V.
 *
 * @tparam Dim The number of dimension for the reflected region
 * @tparam TRegion The region where the I needs to be reflected
 */
template <std::size_t Dim, typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() < 1),
	plsm::SpaceVector<double, Dim>>
getReflectedDispersionForCoefs(const TRegion& clReg)
{
	return clReg.dispersion();
}

/**
 * @brief Computes the sum of cl1RR[i] distance to the mean
 * over the reaction overlap volume.
 *
 * @tparam TRRegion The reflected region
 */
template <typename TRRegion>
KOKKOS_INLINE_FUNCTION
double
computeFirstOrderSum(const typename TRRegion::ScalarType i,
	const TRRegion& cl1RR, const TRRegion& cl2RR, const TRRegion& cl3RR,
	const TRRegion& cl4RR)
{
	double toReturn = 0.0;
	for (auto m : makeIntervalRange(cl4RR[i]))
		for (auto l : makeIntervalRange(cl2RR[i])) {
			toReturn += util::firstOrderSum(
				util::max(cl3RR[i].begin() + m - l, cl1RR[i].begin()),
				util::min(cl3RR[i].end() - 1 + m - l, cl1RR[i].end() - 1),
				static_cast<double>(cl1RR[i].end() - 1 + cl1RR[i].begin()) /
					2.0);
		}
	return toReturn;
}

/**
 * @brief Computes the sum of cl1RR[i] distance to the mean,
 * squared, over the reaction overlap volume.
 *
 * @tparam TRRegion The reflected region
 */
template <typename TRRegion>
KOKKOS_INLINE_FUNCTION
double
computeSecondOrderSum(const typename TRRegion::ScalarType i,
	const TRRegion& cl1RR, const TRRegion& cl2RR, const TRRegion& cl3RR,
	const TRRegion& cl4RR)
{
	double toReturn = 0.0;
	for (auto m : makeIntervalRange(cl4RR[i]))
		for (auto l : makeIntervalRange(cl2RR[i])) {
			toReturn += util::secondOrderSum(
				util::max(cl3RR[i].begin() + m - l, cl1RR[i].begin()),
				util::min(cl3RR[i].end() - 1 + m - l, cl1RR[i].end() - 1),
				static_cast<double>(cl1RR[i].end() - 1 + cl1RR[i].begin()) /
					2.0);
		}
	return toReturn;
}

/**
 * @brief Computes the sum of cl1RR[i] distance to the mean times
 * cl3RR[i] distance to the mean over the reaction overlap volume.
 *
 * @tparam TRRegion The reflected region
 */
template <typename TRRegion>
KOKKOS_INLINE_FUNCTION
double
computeSecondOrderOffsetSum(const typename TRRegion::ScalarType i,
	const TRRegion& cl1RR, const TRRegion& cl2RR, const TRRegion& cl3RR,
	const TRRegion& cl4RR)
{
	double toReturn = 0.0;
	for (auto m : makeIntervalRange(cl4RR[i]))
		for (auto l : makeIntervalRange(cl2RR[i])) {
			toReturn += util::secondOrderOffsetSum(
				util::max(cl3RR[i].begin() + m - l, cl1RR[i].begin()),
				util::min(cl3RR[i].end() - 1 + m - l, cl1RR[i].end() - 1),
				static_cast<double>(cl1RR[i].end() - 1 + cl1RR[i].begin()) /
					2.0,
				static_cast<double>(cl3RR[i].end() - 1 + cl3RR[i].begin()) /
					2.0,
				l - m);
		}
	return toReturn;
}

/**
 * @brief Computes the sum of cl1RR[i] distance to the mean,
 * cubed, over the reaction overlap volume.
 *
 * @tparam TRRegion The reflected region
 */
template <typename TRRegion>
KOKKOS_INLINE_FUNCTION
double
computeThirdOrderSum(const typename TRRegion::ScalarType i,
	const TRRegion& cl1RR, const TRRegion& cl2RR, const TRRegion& cl3RR,
	const TRRegion& cl4RR)
{
	double toReturn = 0.0;
	for (auto m : makeIntervalRange(cl4RR[i]))
		for (auto l : makeIntervalRange(cl2RR[i])) {
			toReturn += pow(l -
								static_cast<double>(
									cl2RR[i].end() - 1 + cl2RR[i].begin()) /
									2.0,
							2.0) *
				util::firstOrderSum(
					util::max(cl3RR[i].begin() + m - l, cl1RR[i].begin()),
					util::min(cl3RR[i].end() - 1 + m - l, cl1RR[i].end() - 1),
					static_cast<double>(cl1RR[i].end() - 1 + cl1RR[i].begin()) /
						2.0);
		}
	return toReturn;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
