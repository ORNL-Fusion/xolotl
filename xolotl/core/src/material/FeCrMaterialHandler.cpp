#include <xolotl/core/material/FeCrMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto feCrMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		FeCrMaterialHandler>({"FeCr"});
}
} // namespace material
} // namespace core
} // namespace xolotl
