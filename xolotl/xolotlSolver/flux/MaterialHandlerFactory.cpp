#include "MaterialHandlerFactory.h"
#include "WFitFluxHandler.h"
#include <iostream>

namespace xolotlSolver {

static std::shared_ptr<IFluxHandler> theFluxHandler;

bool initializeMaterial() {
	bool ret = true;

	// we are to use a tungsten, W, flux handler
	theFluxHandler = std::make_shared<WFitFluxHandler>();

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<IFluxHandler> getMaterialHandler(void) {
	if (!theFluxHandler) {
		// We have not yet been initialized.
		// Issue a warning and use the default tungsten registry.
		std::cerr << "Warning: xolotlSolver material handler requested, but "
				"library has not been initialized.  Defaulting to tungsten "
				"(W) handlers" << std::endl;

		xolotlSolver::initializeMaterial();
	}
	return theFluxHandler;
}

} // end namespace xolotlSolver

