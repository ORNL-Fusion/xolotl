#ifndef TEMPERATUREHANDLERFACTORY_H
#define TEMPERATUREHANDLERFACTORY_H

#include <memory>
#include "ITemperatureHandler.h"
#include <Options.h>

namespace xolotlSolver
{

/**
 * Build the desired type of handler registry.
 * TODO determine if we need to take an enum instead of a bool,
 * if we need to support more than these two types of registries.
 *
 * @param options Options for the program
 * @return True iff the handler registry was created successfully.
 */
bool initializeTempHandler(xolotlCore::Options &options);

/**
 * Access the handler registry.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<ITemperatureHandler> getTemperatureHandler( xolotlCore::Options &options );

}; // end namespace xolotlPerf



#endif /* TEMPERATUREHANDLERFACTORY_H */
