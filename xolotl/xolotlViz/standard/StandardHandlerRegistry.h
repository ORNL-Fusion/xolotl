#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H

#include <iostream>
#include <string>
#include <map>
#include "IVizHandlerRegistry.h"


namespace xolotlViz {

/**
 * Factory for creating plots
 */
class StandardHandlerRegistry : public IVizHandlerRegistry
{
public:

    /**
     * Construct a StandardHandlerRegistry.
     */
    StandardHandlerRegistry();

    /**
     * Clean up a StandardHandlerRegistry.
     */
    virtual ~StandardHandlerRegistry();

    /**
     * Obtain a Plot by name.
     *
     * @param name The name of the Plot.
     * @return A shared pointer to the newly-created Plot.
     */
    virtual std::shared_ptr<IPlot> getPlot(std::string name, PlotType type);

};  //end class StandardHandlerRegistry

}//end namespace xolotlViz

#endif
