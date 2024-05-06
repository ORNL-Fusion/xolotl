#include "../../include/xolotl/core/material/A800H9MeVMaterialHandler.h"

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto a800H9MeVMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		A800H9MeVMaterialHandler>({"800H9MeV"});
}
} // namespace material
} // namespace core
} // namespace xolotl
