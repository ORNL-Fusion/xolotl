#ifndef XCORE_NDPOINT_H
#define XCORE_NDPOINT_H

#include <array>

namespace xolotl {
namespace util {

/**
 * Nice name for coordinates of a 3D point.
 */
template<uint32_t N>
using Point = std::array<double, N>;

} /* end namespace util */
} /* end namespace xolotl */

#endif // XCORE_NDPOINT_H
