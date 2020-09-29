#include <xolotl/factory/viz/VizHandlerFactory.h>
#include <xolotl/viz/dummy/DummyHandler.h>

namespace xolotl
{
namespace viz
{
namespace dummy
{
namespace detail
{
auto dummyHandlerRegistrations =
	::xolotl::factory::viz::VizHandlerFactory::RegistrationCollection<
		DummyHandler>({"dummy"});
}

DummyHandler::DummyHandler(const options::IOptions&)
{
}

DummyHandler::~DummyHandler() = default;

std::shared_ptr<IPlot>
DummyHandler::getPlot(const std::string& name, PlotType)
{
	return std::make_shared<DummyPlot>(name);
}

} // namespace dummy
} // namespace viz
} // namespace xolotl
