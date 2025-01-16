#pragma once

#include <Kokkos_Atomic.hpp>

namespace xolotl
{
namespace util
{
template <typename T>
KOKKOS_INLINE_FUNCTION
bool
atomicCompareExchangeStrong(T* ptr, T expected, T desired)
{
    return expected == Kokkos::atomic_compare_exchange(ptr, expected, desired);
}
}
}
