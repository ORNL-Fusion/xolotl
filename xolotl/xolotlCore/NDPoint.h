#ifndef XCORE_NDPOINT_H
#define XCORE_NDPOINT_H

#include <array>

namespace xolotlCore {

/**
 * Nice name for coordinates of a 3D point.
 */
template<uint32_t N>
using NDPoint = std::array<double, N>;

} // namespace xolotlCore

#endif // XCORE_NDPOINT_H
