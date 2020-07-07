#pragma once

#include <xolotl/factory/viz/VizHandlerRegistryFactory.h>
#include <xolotl/options/Options.h>
#include <xolotl/viz/IVizHandlerRegistry.h>

namespace xolotl
{
namespace viz
{
class VizHandlerRegistry : public IVizHandlerRegistry
{
public:
	VizHandlerRegistry(const options::Options& options);

	std::shared_ptr<IPlot>
	getPlot(const std::string& name, PlotType type);

	static IVizHandlerRegistry*
	get()
	{
		return staticVizHandlerRegistry.get();
	}

private:
	static std::shared_ptr<IVizHandlerRegistry> staticVizHandlerRegistry;
};
} // namespace viz
} // namespace xolotl
