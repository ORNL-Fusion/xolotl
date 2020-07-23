#include <xolotl/core/material/FuelMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto fuelMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		FuelMaterialHandler>({"Fuel"});
}
} // namespace material
} // namespace core
} // namespace xolotl
