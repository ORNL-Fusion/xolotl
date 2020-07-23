#include <xolotl/core/material/AlloyMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto alloyMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		AlloyMaterialHandler>({"800H"});
}
} // namespace material
} // namespace core
} // namespace xolotl
