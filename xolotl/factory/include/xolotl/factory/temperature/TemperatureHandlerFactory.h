#ifndef TEMPERATUREHANDLERFACTORY_H
#define TEMPERATUREHANDLERFACTORY_H

#include <memory>

#include <xolotl/core/temperature/ITemperatureHandler.h>
#include <xolotl/options/Options.h>

namespace xolotl
{
namespace factory
{
namespace temperature
{
/**
 * Build the desired type of handler registry.
 *
 * @param options Options for the program
 * @return True iff the handler registry was created successfully.
 */
bool
initializeTempHandler(const options::Options& opts);

/**
 * Access the handler registry.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<core::temperature::ITemperatureHandler>
getTemperatureHandler();

} // end namespace temperature
} // end namespace factory
} // end namespace xolotl

#endif /* TEMPERATUREHANDLERFACTORY_H */
