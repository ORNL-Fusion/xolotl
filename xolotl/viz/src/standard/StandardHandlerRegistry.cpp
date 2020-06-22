#include <xolotl/viz/LabelProvider.h>
#include <xolotl/viz/dataprovider/DataProvider.h>
#include <xolotl/viz/standard/StandardHandlerRegistry.h>
#include <xolotl/viz/standard/plot/Plot.h>
#include <xolotl/viz/standard/plot/ScatterPlot.h>
#include <xolotl/viz/standard/plot/SeriesPlot.h>
#include <xolotl/viz/standard/plot/SurfacePlot.h>
#include <xolotl/viz/standard/plot/VideoPlot.h>

namespace xolotl
{
namespace viz
{
namespace standard
{
StandardHandlerRegistry::StandardHandlerRegistry()
{
}

StandardHandlerRegistry::~StandardHandlerRegistry()
{
}

std::shared_ptr<IPlot>
StandardHandlerRegistry::getPlot(const std::string& name, PlotType type)
{
	switch (type) {
	case PlotType::SCATTER:
		return std::make_shared<plot::ScatterPlot>(name);
	case PlotType::SERIES:
		return std::make_shared<plot::SeriesPlot>(name);
	case PlotType::SURFACE:
		return std::make_shared<plot::SurfacePlot>(name);
	case PlotType::VIDEO:
		return std::make_shared<plot::VideoPlot>(name);
	default:
		return std::make_shared<plot::Plot>(name);
	}
}

} // namespace standard
} // namespace viz
} // namespace xolotl
