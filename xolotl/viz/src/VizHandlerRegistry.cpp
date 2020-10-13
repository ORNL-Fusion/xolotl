#include <stdexcept>

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

const std::shared_ptr<IVizHandler>&
VizHandlerRegistry::get()
{
	if (!vizHandler) {
		throw std::runtime_error(
			"Request for viz handler before it has been registered");
	}
	return vizHandler;
}
} // namespace viz
} // namespace xolotl
