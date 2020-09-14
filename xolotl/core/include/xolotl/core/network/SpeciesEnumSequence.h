#pragma once

#include <string>
#include <type_traits>

#include <xolotl/core/network/EnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TSpeciesEnum>
std::string
getSpeciesLabel(TSpeciesEnum species);

template <typename TSpeciesEnum>
struct NumberOfInterstitialSpecies : std::integral_constant<std::size_t, 0>
{
};

template <typename TSpeciesEnum>
struct NumberOfVacancySpecies : std::integral_constant<std::size_t, 0>
{
};

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr std::size_t
numberOfInterstitialSpecies() noexcept
{
	return NumberOfInterstitialSpecies<TSpeciesEnum>::value;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr std::size_t
numberOfVacancySpecies() noexcept
{
	return NumberOfVacancySpecies<TSpeciesEnum>::value;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr bool
isVacancy(TSpeciesEnum val,
	std::enable_if_t<(numberOfInterstitialSpecies<TSpeciesEnum>() > 1), int> =
		0) noexcept
{
	return val == TSpeciesEnum::V || val == TSpeciesEnum::Void ||
		val == TSpeciesEnum::Faulted;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr bool
isVacancy(TSpeciesEnum val,
	std::enable_if_t<numberOfInterstitialSpecies<TSpeciesEnum>() == 1, int> =
		0) noexcept
{
	return val == TSpeciesEnum::V;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr bool
isVacancy(TSpeciesEnum val,
	std::enable_if_t<numberOfInterstitialSpecies<TSpeciesEnum>() == 0, int> =
		0) noexcept
{
	return false;
}

template <typename TSpeciesEnum, std::size_t N>
KOKKOS_INLINE_FUNCTION
constexpr bool
isVacancy(EnumSequence<TSpeciesEnum, N> seq) noexcept
{
	return isVacancy(seq.value);
}

template <typename TSpeciesEnum, std::size_t N>
class SpeciesEnumSequence : public EnumSequence<TSpeciesEnum, N>
{
public:
	using Sequence = EnumSequence<TSpeciesEnum, N>;

	using Sequence::Sequence;

	KOKKOS_INLINE_FUNCTION
	static constexpr std::size_t
	sizeNoI() noexcept
	{
		return numberOfVacancySpecies<TSpeciesEnum>() > 0 ?
			N - numberOfInterstitialSpecies<TSpeciesEnum>() -
				numberOfVacancySpecies<TSpeciesEnum>() + 1 :
			N;
	}

	KOKKOS_INLINE_FUNCTION
	static constexpr Sequence
	lastNoI() noexcept
	{
		return numberOfInterstitialSpecies<TSpeciesEnum>() > 0 ?
			Sequence::last() - numberOfInterstitialSpecies<TSpeciesEnum>() -
				numberOfVacancySpecies<TSpeciesEnum>() + 1 :
			Sequence::last();
	}
};

template <typename TSpeciesEnum, std::size_t N>
struct SpeciesForGrouping
{
	using Sequence = SpeciesEnumSequence<TSpeciesEnum, N>;
	static constexpr auto first = Sequence::first();
	static constexpr auto last = Sequence::lastNoI();

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<TSpeciesEnum>
	mapToMomentId(EnumSequence<TSpeciesEnum, N> value)
	{
		return value();
	}
};
} // namespace network
} // namespace core
} // namespace xolotl
