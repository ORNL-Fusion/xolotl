#include <xolotl/util/MPIUtils.h>
#include <xolotl/viz/VizHandlerRegistry.h>
#include <xolotl/viz/config.h>

namespace xolotl
{
namespace viz
{
std::shared_ptr<IVizHandler> VizHandlerRegistry::vizHandler;

void
VizHandlerRegistry::set(const std::shared_ptr<IVizHandler>& handler)
{
	vizHandler = handler;
}
} // namespace viz
} // namespace xolotl
