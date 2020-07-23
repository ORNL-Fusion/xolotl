#include <xolotl/core/material/PulsedMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto pulsedMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		PulsedMaterialHandler>({"Pulsed"});
}
} // namespace material
} // namespace core
} // namespace xolotl
