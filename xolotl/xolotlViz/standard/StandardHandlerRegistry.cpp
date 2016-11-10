#include "StandardHandlerRegistry.h"
#include <Plot.h>
#include <ScatterPlot.h>
#include <SeriesPlot.h>
#include <SurfacePlot.h>
#include <VideoPlot.h>
#include <DataProvider.h>
#include <LabelProvider.h>

namespace xolotlViz {

StandardHandlerRegistry::StandardHandlerRegistry() {
}

StandardHandlerRegistry::~StandardHandlerRegistry() {
}

std::shared_ptr<IPlot> StandardHandlerRegistry::getPlot(const std::string& name,
		PlotType type) {
	switch (type) {
	case PlotType::SCATTER:
		return std::make_shared<ScatterPlot>(name);
	case PlotType::SERIES:
		return std::make_shared<SeriesPlot>(name);
	case PlotType::SURFACE:
		return std::make_shared<SurfacePlot>(name);
	case PlotType::VIDEO:
		return std::make_shared<VideoPlot>(name);
	default:
		return std::make_shared<Plot>(name);
	}
}

}    //end namespace xolotlViz

