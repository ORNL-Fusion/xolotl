#include "MaterialHandlerFactory.h"
#include "WFitFluxHandler.h"
#include <FluxHandler.h>
#include <iostream>

namespace xolotlSolver {

static std::shared_ptr<IFluxHandler> theFluxHandler;

bool initializeMaterial( xolotlCore::Options &options ) {
	bool ret = true;

	// we are to use a tungsten, W, flux handler
	theFluxHandler = std::make_shared<WFitFluxHandler>();

	// If the Helium fluence option is present, set the value
	if( options.useMaxHeliumFluence() )
	{
		theFluxHandler->setMaxHeFluence( options.getMaxHeliumFluence() );
	}

	// If the helium flux option is present, set the value
	if( options.useHeliumFlux() )
	{
		theFluxHandler->setHeFlux( options.getHeliumFlux() );
	}

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<IFluxHandler> getMaterialHandler(void) {
	if (!theFluxHandler) {
		// We have not yet been initialized.
		// Issue a warning.
		std::cerr << "Warning: xolotlSolver material handler requested, but "
				"library has not been initialized." << std::endl;
	}

	return theFluxHandler;
}

} // end namespace xolotlSolver

