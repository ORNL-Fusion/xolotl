#pragma once

// Includes
#include <memory>
#include <string>

#include <xolotl/viz/IPlot.h>
#include <xolotl/viz/PlotType.h>

namespace xolotl
{
namespace viz
{
/**
 * Interface for visualization handlers
 */
class IVizHandler
{
public:
	/**
	 * The destructor
	 */
	virtual ~IVizHandler()
	{
	}

	/**
	 * Obtain a Plot by type.
	 *
	 * @param type The type of plot to return.
	 * @return A shared pointer to the newly-created Plot.
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(PlotType type) = 0;

}; // end class IVizHandler

} /* namespace viz */
} /* namespace xolotl */
