#ifndef DUMMYHANDLERREGISTRY_H
#define DUMMYHANDLERREGISTRY_H

#include "IVizHandlerRegistry.h"
#include <DummyPlot.h>

namespace xolotlViz {

/**
 * Factory for creating plots
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
    virtual std::shared_ptr<IPlot> getPlot(std::string name, PlotType type);

};  //end class DummyHandlerRegistry

}//end namespace xolotlViz

#endif
