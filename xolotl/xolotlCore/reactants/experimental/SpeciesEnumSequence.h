#pragma once

#include <type_traits>

#include <experimental/EnumSequence.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
struct HasInterstitial : std::false_type
{
};

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr bool
hasInterstitial() noexcept
{
    return HasInterstitial<TSpeciesEnum>::value;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr bool
isVacancy(TSpeciesEnum val,
    std::enable_if_t<hasInterstitial<TSpeciesEnum>(), int> = 0) noexcept
{
    return val == TSpeciesEnum::V;
}

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr bool
isVacancy(TSpeciesEnum val,
    std::enable_if_t<!hasInterstitial<TSpeciesEnum>(), int> = 0) noexcept
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
    static
    constexpr std::size_t
    sizeNoI() noexcept
    {
        return hasInterstitial<TSpeciesEnum>() ? N - 1 : N;
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr Sequence
    lastNoI() noexcept
    {
        return hasInterstitial<TSpeciesEnum>() ?
            Sequence::last() - 1 : Sequence::last();
    }
};
}
}
