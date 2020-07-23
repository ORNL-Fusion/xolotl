#include <xolotl/core/material/W211MaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto w211MaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		W211MaterialHandler>({"W211"});
}
} // namespace material
} // namespace core
} // namespace xolotl
