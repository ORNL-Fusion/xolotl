#pragma once

#include <Kokkos_View.hpp>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
using DefaultMemorySpace = typename Kokkos::View<int*>::traits::memory_space;
}
}
}
