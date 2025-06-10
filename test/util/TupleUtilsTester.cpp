#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>

#include <xolotl/util/TupleUtils.h>

BOOST_AUTO_TEST_SUITE(TupleUtils_testSuite)

using namespace xolotl::util;
struct A
{
};
struct B
{
};
struct C
{
};
template <typename... T>
struct WrapAll
{
};

BOOST_AUTO_TEST_CASE(TupleUtils)
{
	using T0 = std::tuple<>;

	// TuplePushFront
	using T1 = TuplePushFront<T0, A>;
	static_assert(std::is_same_v<T1, std::tuple<A>>);
	using T2 = TuplePushFront<T1, B>;
	static_assert(std::is_same_v<T2, std::tuple<B, A>>);
	using T3 = TuplePushFront<T2, C>;
	static_assert(std::is_same_v<T3, std::tuple<C, B, A>>);

	// TuplePopFront
	using T4 = TuplePopFront<T3>;
	static_assert(std::is_same_v<T4, T2>);
	using T5 = TuplePopFront<T4>;
	static_assert(std::is_same_v<T5, T1>);
	using T6 = TuplePopFront<T5>;
	static_assert(std::is_same_v<T6, T0>);
	using T60 = TuplePopFront<T6>;
	static_assert(std::is_same_v<T60, T0>);

	// TuplePushBack
	using T7 = TuplePushBack<T6, A>;
	static_assert(std::is_same_v<T7, std::tuple<A>>);
	using T8 = TuplePushBack<T7, B>;
	static_assert(std::is_same_v<T8, std::tuple<A, B>>);
	using T9 = TuplePushBack<T8, C>;
	static_assert(std::is_same_v<T9, std::tuple<A, B, C>>);

	// TuplePopBack
	using T10 = TuplePopBack<T9>;
	static_assert(std::is_same_v<T10, T8>);
	using T11 = TuplePopBack<T10>;
	static_assert(std::is_same_v<T11, T7>);
	using T12 = TuplePopBack<T11>;
	static_assert(std::is_same_v<T12, T0>);
	using T120 = TuplePopBack<T12>;
	static_assert(std::is_same_v<T120, T0>);

	// TupleReverse
	using T13 = TupleReverse<T9>;
	static_assert(std::is_same_v<T13, T3>);
	using T14 = TupleReverse<T13>;
	static_assert(std::is_same_v<T14, T9>);

	// TupleApplyAll
	using T15 = TupleApplyAll<WrapAll, T3>;
	static_assert(std::is_same_v<T15, WrapAll<C, B, A>>);
}

BOOST_AUTO_TEST_SUITE_END()
