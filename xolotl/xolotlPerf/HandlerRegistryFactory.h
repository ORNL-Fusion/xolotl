#ifndef HANDLERREGISTRYFACTORY_H
#define HANDLERREGISTRYFACTORY_H

#include <memory>
#include <vector>
#include "IHandlerRegistry.h"
#include "HardwareQuantities.h"


namespace xolotlPerf
{

// Build the desired type of handler registry.
// 
// TODO determine if we need to take an enum instead of a bool,
// if we need to support more than these two types of registries.
//
// @param useStdRegistry Whether to use the "standard" handlers or dummy handlers.
// @param hwQuantities Collection of hardware quantities to monitor along with our timings.  (Default: none)
// @return True iff the handler registry was created successfully.
bool initialize( bool useStdRegistry,
                    std::vector<HardwareQuantities> hwQuantities = std::vector<HardwareQuantities>() );

// Access the handler registry.
// @return The handler registry object.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry( void );

}; // end namespace xolotlPerf

#endif // HANDLERREGISTRYFACTORY_H
