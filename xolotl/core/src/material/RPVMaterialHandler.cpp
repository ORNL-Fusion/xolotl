#include <xolotl/core/material/RPVMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto rpvMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		RPVMaterialHandler>({"RPV"});
}
} // namespace material
} // namespace core
} // namespace xolotl
