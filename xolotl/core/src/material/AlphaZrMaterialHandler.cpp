#include <xolotl/core/material/AlphaZrMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto zrMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		AlphaZrMaterialHandler>({"AlphaZr"});
}
} // namespace material
} // namespace core
} // namespace xolotl
