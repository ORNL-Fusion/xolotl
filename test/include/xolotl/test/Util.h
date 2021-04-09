#pragma once

#include <cmath>

namespace xolotl::test
{
inline constexpr double defaultZeroTol = 1.0e-12;

} // namespace xolotl::test

#define XOLOTL_REQUIRE_CLOSE_ZT(a, b, relTol, zeroTol) \
	{ \
		if (std::fabs(a) < zeroTol || std::fabs(b) < zeroTol) { \
			BOOST_REQUIRE_SMALL(std::fabs(a - b), zeroTol); \
		} \
		else { \
			BOOST_REQUIRE_CLOSE(a, b, relTol); \
		} \
	}

#define XOLOTL_REQUIRE_CLOSE(a, b, relTol) \
	XOLOTL_REQUIRE_CLOSE_ZT(a, b, relTol, ::xolotl::test::defaultZeroTol)
