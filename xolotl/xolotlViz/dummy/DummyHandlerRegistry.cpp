#include "DummyHandlerRegistry.h"

namespace xolotlViz {

DummyHandlerRegistry::DummyHandlerRegistry() {
}

DummyHandlerRegistry::~DummyHandlerRegistry() {
}

std::shared_ptr<IPlot> DummyHandlerRegistry::getPlot(const std::string& name, PlotType) {
	return std::make_shared <DummyPlot> (name);
}

}    //end namespace xolotlViz

