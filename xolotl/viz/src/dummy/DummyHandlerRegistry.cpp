#include <xolotl/viz/dummy/DummyHandlerRegistry.h>

namespace xolotl
{
namespace viz
{
namespace dummy
{
DummyHandlerRegistry::DummyHandlerRegistry()
{
}

DummyHandlerRegistry::~DummyHandlerRegistry()
{
}

std::shared_ptr<IPlot>
DummyHandlerRegistry::getPlot(const std::string& name, PlotType)
{
	return std::make_shared<DummyPlot>(name);
}

} // namespace dummy
} // namespace viz
} // namespace xolotl
