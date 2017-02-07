#ifndef IVIZHANDLERREGISTRY_H
#define IVIZHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "IPlot.h"
#include "PlotType.h"


namespace xolotlViz {

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

} //end namespace xolotlViz

#endif
