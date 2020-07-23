#include <xolotl/core/material/FeMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto feMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		FeMaterialHandler>({"Fe"});
}
} // namespace material
} // namespace core
} // namespace xolotl
