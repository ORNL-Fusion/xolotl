#pragma once

// Includes
#include <memory>
#include <string>
#include <vector>

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
	 * This operation returns the IPlot specified by the parameter.
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(const std::string& name, PlotType type) = 0;

}; // end class IVizHandler

} /* namespace viz */
} /* namespace xolotl */
