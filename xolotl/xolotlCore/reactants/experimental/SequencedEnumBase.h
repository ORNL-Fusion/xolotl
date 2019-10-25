#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TDerived, std::size_t NumElements>
class SequencedEnumBase
{
    // template <std::size_t... Is>
    // static constexpr auto
    // makeIndexTuple(std::index_sequence<Is...>)
    // {
    //     return std::make_tuple(Is...);
    // }

public:
    // static constexpr auto indices =
    //     makeIndexTuple(std::make_index_sequence<NumElements>{});

    constexpr SequencedEnumBase() : value{-1} { }
    constexpr SequencedEnumBase(int i) : value{i} { }
    constexpr SequencedEnumBase(const SequencedEnumBase&) = default;

    constexpr
    operator TDerived() const noexcept
    {
        return value;
    }

    // explicit
    constexpr
    operator std::size_t() const noexcept
    {
        return value;
    }

    static constexpr
    TDerived
    first() noexcept
    {
        return TDerived{0};
    }

    static constexpr
    TDerived
    last() noexcept
    {
        return TDerived{NumElements - 1};
    }

    static constexpr
    TDerived
    invalid() noexcept
    {
        return TDerived{};
    }

    //TODO: Should this wrap?
    SequencedEnumBase&
    operator++()
    {
        ++value;
        return *this;
    }

    int value;
};


template <typename T, std::size_t N>
inline
constexpr bool
operator==(SequencedEnumBase<T, N> a, SequencedEnumBase<T, N> b)
{
    return a.value == b.value;
}


template <typename T, std::size_t N>
inline
constexpr bool
operator!=(SequencedEnumBase<T, N> a, SequencedEnumBase<T, N> b)
{
    return !(a == b);
}


template <typename T, std::size_t N>
inline
constexpr bool
operator<(SequencedEnumBase<T, N> a, SequencedEnumBase<T, N> b)
{
    return a.value < b.value;
}


template <typename T, std::size_t N>
inline
constexpr bool
operator<=(SequencedEnumBase<T, N> a, SequencedEnumBase<T, N> b)
{
    return a.value <= b.value;
}


template <typename T, std::size_t N>
inline
constexpr bool
operator>(SequencedEnumBase<T, N> a, SequencedEnumBase<T, N> b)
{
    return a.value > b.value;
}


template <typename T, std::size_t N>
inline
constexpr bool
operator>=(SequencedEnumBase<T, N> a, SequencedEnumBase<T, N> b)
{
    return a.value >= b.value;
}
}
}
