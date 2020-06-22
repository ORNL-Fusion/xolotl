#ifndef DUMMYHANDLERREGISTRY_H
#define DUMMYHANDLERREGISTRY_H

#include <xolotl/viz/IVizHandlerRegistry.h>
#include <xolotl/viz/dummy/DummyPlot.h>

namespace xolotl
{
namespace viz
{
namespace dummy
{
/**
 * Factory for creating plots that are dummies, meaning that they have the
 * same structure as IPlot but don't actually do anything. This is so that
 * the code can be written to use the visualization infrastructure without
 * regard to whether visualization is active or disabled.
 */
class DummyHandlerRegistry : public IVizHandlerRegistry
{
public:
	/**
	 * Construct a DummyHandlerRegistry.
	 */
	DummyHandlerRegistry();

	/**
	 * Clean up a DummyHandlerRegistry.
	 */
	virtual ~DummyHandlerRegistry();

	/**
	 * Obtain a Plot by name.
	 *
	 * @param name The name of the Plot.
	 * @return A shared pointer to the newly-created Plot.
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(const std::string& name, PlotType type);
};
// end class DummyHandlerRegistry

} // end namespace dummy
} // end namespace viz
} // end namespace xolotl

#endif
