#include <xolotl/core/material/W111MaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto w111MaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		W111MaterialHandler>({"W111"});
}
} // namespace material
} // namespace core
} // namespace xolotl
