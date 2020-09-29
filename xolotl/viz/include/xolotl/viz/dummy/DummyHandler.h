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

	/**
	 * Clean up a DummyHandler.
	 */
	virtual ~DummyHandler();

	/**
	 * Obtain a Plot by name.
	 *
	 * @param name The name of the Plot.
	 * @return A shared pointer to the newly-created Plot.
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(const std::string& name, PlotType type);
};
// end class DummyHandler

} // end namespace dummy
} // end namespace viz
} // end namespace xolotl
