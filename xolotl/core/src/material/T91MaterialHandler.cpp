#include <xolotl/core/material/T91MaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto t91MaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		T91MaterialHandler>({"T91"});
}
} // namespace material
} // namespace core
} // namespace xolotl
