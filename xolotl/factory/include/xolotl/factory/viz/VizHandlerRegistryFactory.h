#ifndef VIZHANDLERREGISTRYFACTORY_H
#define VIZHANDLERREGISTRYFACTORY_H

#include <memory>
#include <vector>
#include <xolotl/viz/IVizHandlerRegistry.h>


namespace xolotl
{
namespace factory
{
namespace viz
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
std::shared_ptr<::xolotl::viz::IVizHandlerRegistry> getVizHandlerRegistry();

} // end namespace viz
} // end namespace factory
} // end namespace xolotl

#endif // VIZHANDLERREGISTRYFACTORY_H
