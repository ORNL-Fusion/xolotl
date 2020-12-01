#pragma once

#include <xolotl/viz/IVizHandler.h>
#include <xolotl/viz/dummy/DummyPlot.h>

namespace xolotl
{
namespace viz
{
namespace dummy
{
/**
 * Handler for creating plots that are dummies, meaning that they have the
 * same structure as IPlot but don't actually do anything. This is so that
 * the code can be written to use the visualization infrastructure without
 * regard to whether visualization is active or disabled.
 */
class DummyHandler : public IVizHandler
{
public:
	DummyHandler() = delete;

	DummyHandler(const options::IOptions& options);

	virtual ~DummyHandler();

	/**
	 * \see IVizHandler.h
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(PlotType type);
};
// end class DummyHandler

} // end namespace dummy
} // end namespace viz
} // end namespace xolotl
