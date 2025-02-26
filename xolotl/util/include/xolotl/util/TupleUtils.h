#pragma once

#include <tuple>

namespace xolotl
{
namespace util
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
} // namespace detail

template <typename TTuple, typename T>
using TuplePushFront = typename detail::TuplePushFrontHelper<T, TTuple>::Type;

namespace detail
{
template <typename TTuple>
struct TuplePopFrontHelper;

template <>
struct TuplePopFrontHelper<std::tuple<>>
{
	using Type = std::tuple<>;
};

template <typename T, typename... Ts>
struct TuplePopFrontHelper<std::tuple<T, Ts...>>
{
	using Type = std::tuple<Ts...>;
};
} // namespace detail

template <typename TTuple>
using TuplePopFront = typename detail::TuplePopFrontHelper<TTuple>::Type;

namespace detail
{
template <typename TTuple, typename T>
struct TuplePushBackHelper;

template <typename... Ts, typename T>
struct TuplePushBackHelper<std::tuple<Ts...>, T>
{
	using Type = std::tuple<Ts..., T>;
};
} // namespace detail

template <typename TTuple, typename T>
using TuplePushBack = typename detail::TuplePushBackHelper<TTuple, T>::Type;

namespace detail
{
template <typename TTuple>
struct TuplePopBackHelper;

template <>
struct TuplePopBackHelper<std::tuple<>>
{
	using Type = std::tuple<>;
};

template <typename T>
struct TuplePopBackHelper<std::tuple<T>>
{
	using Type = std::tuple<>;
};

template <typename T, typename... Ts>
struct TuplePopBackHelper<std::tuple<T, Ts...>>
{
	using Type =
		TuplePushFront<typename TuplePopBackHelper<std::tuple<Ts...>>::Type, T>;
};
} // namespace detail

template <typename TTuple>
using TuplePopBack = typename detail::TuplePopBackHelper<TTuple>::Type;

namespace detail
{
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
} // namespace detail

template <typename TTuple>
using TupleReverse = typename detail::TupleReverseHelper<TTuple>::Type;

namespace detail
{
template <template <typename...> typename Tpl, typename TTuple>
struct TupleApplyAllHelper;

template <template <typename...> typename Tpl, typename... Ts>
struct TupleApplyAllHelper<Tpl, std::tuple<Ts...>>
{
	using Type = Tpl<Ts...>;
};
} // namespace detail

template <template <typename...> typename Tpl, typename TTuple>
using TupleApplyAll = typename detail::TupleApplyAllHelper<Tpl, TTuple>::Type;
} // namespace util
} // namespace xolotl
