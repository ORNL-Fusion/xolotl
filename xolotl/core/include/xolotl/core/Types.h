#include <Kokkos_Array.hpp>
#include <Kokkos_View.hpp>

#include <xolotl/config.h>

namespace xolotl
{
namespace core
{
using StencilConcArray = Kokkos::Array<Kokkos::View<const double*>,
	KOKKOS_INVALID_INDEX, Kokkos::Array<>::contiguous>;
}
} // namespace xolotl
