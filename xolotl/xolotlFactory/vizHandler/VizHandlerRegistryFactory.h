#ifndef VIZHANDLERREGISTRYFACTORY_H
#define VIZHANDLERREGISTRYFACTORY_H

#include <memory>
#include <vector>
#include <IVizHandlerRegistry.h>


namespace xolotlFactory
{

/**
 * Build the desired type of handler registry.
 * @return True iff the handler registry was created successfully.
 */
bool initializeVizHandler(bool useStdRegistry);

/**
 * Access the handler registry.
 * @return The handler registry object.
 */
std::shared_ptr<xolotlViz::IVizHandlerRegistry> getVizHandlerRegistry();

}; // end namespace xolotlFactory

#endif // VIZHANDLERREGISTRYFACTORY_H
