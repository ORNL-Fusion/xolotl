#ifndef IVIZHANDLERREGISTRY_H
#define IVIZHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <xolotl/viz/IPlot.h>
#include <xolotl/viz/PlotType.h>


namespace xolotl {
namespace viz {

/**
 * Factory for building plots, dataprovider, and labelprovider.
 */
class IVizHandlerRegistry {

public:

	/**
	 * The destructor
	 */
	virtual ~IVizHandlerRegistry(){}

	/**
	 * This operation returns the IPlot specified by the parameter.
	 */
	virtual std::shared_ptr<IPlot> getPlot(const std::string& name, PlotType type) = 0;

}; //end class IVizHandlerRegistry

} /* namespace viz */
} /* namespace xolotl */

#endif
