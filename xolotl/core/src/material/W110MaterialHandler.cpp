#include <xolotl/core/material/W110MaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto w110MaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		W110MaterialHandler>({"W110"});
}
} // namespace material
} // namespace core
} // namespace xolotl
