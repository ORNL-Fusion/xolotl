#include "../../include/xolotl/core/material/A800H5MeVMaterialHandler.h"

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto a800H5MeVMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		A800H5MeVMaterialHandler>({"800H5MeV"});
}
} // namespace material
} // namespace core
} // namespace xolotl
