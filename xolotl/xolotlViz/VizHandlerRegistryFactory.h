#ifndef VIZHANDLERREGISTRYFACTORY_H
#define VIZHANDLERREGISTRYFACTORY_H

#include <memory>
#include <vector>
#include "IVizHandlerRegistry.h"


namespace xolotlViz
{

/**
 * Build the desired type of handler registry.
 * @return True iff the handler registry was created successfully.
 */
bool initialize(bool useStdRegistry);

/**
 * Access the handler registry.
 * @return The handler registry object.
 */
std::shared_ptr<IVizHandlerRegistry> getVizHandlerRegistry();

}; // end namespace xolotlViz

#endif // VIZHANDLERREGISTRYFACTORY_H
