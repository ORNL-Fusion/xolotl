#ifndef TEMPERATUREHANDLERFACTORY_H
#define TEMPERATUREHANDLERFACTORY_H

#include <memory>
#include "ITemperatureHandler.h"
#include <XolotlOptions.h>

namespace xolotlSolver
{

/**
 * Build the desired type of handler registry.
 * TODO determine if we need to take an enum instead of a bool,
 * if we need to support more than these two types of registries.
 *
 * @param useConstRegistry Whether to use the constant temperature handlers.
 * @param useTempProfileRegistry Whether to use the temperature profile handlers.
 * @param options Options for the program
 * @return True iff the handler registry was created successfully.
 */
bool initializeTempHandler( bool useConstTempRegistry, bool useTempProfileRegistry,
		xolotlCore::XolotlOptions &options);

/**
 * Access the handler registry.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<ITemperatureHandler> getTemperatureHandler( xolotlCore::XolotlOptions &options );

}; // end namespace xolotlPerf



#endif /* TEMPERATUREHANDLERFACTORY_H */
