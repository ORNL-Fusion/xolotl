#pragma once

#include <Kokkos_View.hpp>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
using DefaultMemorySpace = typename Kokkos::View<int*>::traits::memory_space;
}
}
}
}
