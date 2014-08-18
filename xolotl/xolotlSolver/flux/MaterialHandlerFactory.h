#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include "IFluxHandler.h"
#include <Options.h>

namespace xolotlSolver
{

/**
 * Build the desired type of handler registry.
 * TODO determine if we need to take an enum instead of a bool,
 * if we need to support more than these two types of registries.
 *
 * @param usesMaxHeFluence	Whether the maximum Helium fluence option was specified
 * @param options	Xolotl options set from the parameter file
 * @return True iff the handler registry was created successfully.
 */
bool initializeMaterial( xolotlCore::Options &options );

//bool initializeMaterial( bool useWRegistry );

/**
 * Access the handler registry.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<IFluxHandler> getMaterialHandler( void );

}; // end namespace xolotlPerf

#endif // MATERIALHANDLERFACTORY_H
