#include <xolotl/core/material/W100MaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto w100MaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		W100MaterialHandler>({"W100"});
}
} // namespace material
} // namespace core
} // namespace xolotl
