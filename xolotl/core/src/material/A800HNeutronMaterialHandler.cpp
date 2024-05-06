#include "../../include/xolotl/core/material/A800HNeutronMaterialHandler.h"

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto a800HNeutronMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		A800HNeutronMaterialHandler>({"800HNeutron"});
}
} // namespace material
} // namespace core
} // namespace xolotl
