#pragma once

#include <tuple>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename T, typename TTuple>
struct TuplePushFrontHelper;

template <typename T, typename... Ts>
struct TuplePushFrontHelper<T, std::tuple<Ts...>>
{
    using Type = std::tuple<T, Ts...>;
};

template <typename T, typename TTuple>
using TuplePushFront = typename TuplePushFrontHelper<T, TTuple>::Type;

template <typename TTuple>
struct TuplePopFrontHelper;

template <typename T, typename... Ts>
struct TuplePopFrontHelper<std::tuple<T, Ts...>>
{
    using Type = std::tuple<Ts...>;
};

template <typename TTuple>
using TuplePopFront = typename TuplePopFrontHelper<TTuple>::Type;

template <typename TTuple, typename T>
struct TuplePushBackHelper;

template <typename... Ts, typename T>
struct TuplePushBackHelper<std::tuple<Ts...>, T>
{
    using Type = std::tuple<Ts..., T>;
};

template <typename TTuple, typename T>
using TuplePushBack = typename TuplePushBackHelper<TTuple, T>::Type;

template <typename TTuple>
struct TuplePopBackHelper;

template <typename T>
struct TuplePopBackHelper<std::tuple<T>>
{
    using Type = std::tuple<>;
};

template <typename T, typename... Ts>
struct TuplePopBackHelper<std::tuple<T, Ts...>>
{
    using Type = TuplePushFront<T,
        typename TuplePopBackHelper<std::tuple<Ts...>>::Type>;
};

template <typename TTuple>
using TuplePopBack = typename TuplePopBackHelper<TTuple>::Type;

template <typename TTuple>
struct TupleReverseHelper;

template <>
struct TupleReverseHelper<std::tuple<>>
{
    using Type = std::tuple<>;
};

template <typename T, typename... Ts>
struct TupleReverseHelper<std::tuple<T, Ts...>>
{
    using Type =
        TuplePushBack<typename TupleReverseHelper<std::tuple<Ts...>>::Type, T>;
};

template <typename TTuple>
using TupleReverse = typename TupleReverseHelper<TTuple>::Type;
}
}
}
}
