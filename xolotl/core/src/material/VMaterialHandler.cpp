#include <xolotl/core/material/VMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto vMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		VMaterialHandler>({"V"});
}
} // namespace material
} // namespace core
} // namespace xolotl
