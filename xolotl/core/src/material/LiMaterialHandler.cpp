#include <xolotl/core/material/LiMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto liMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		LiMaterialHandler>({"Li"});
}
} // namespace material
} // namespace core
} // namespace xolotl
