#include "MaterialHandlerFactory.h"
#include "W100FitFluxHandler.h"
#include "W111FitFluxHandler.h"
#include "FeFitFluxHandler.h"
#include <FluxHandler.h>
#include <iostream>

namespace xolotlSolver {

static std::shared_ptr<IFluxHandler> theFluxHandler;

bool initializeMaterial(xolotlCore::Options &options) {
	bool ret = true;

	if (options.useMaterial()) {
		std::string materialName = options.getMaterial();

		// W (100) case
		if (materialName == "W100") {
			theFluxHandler = std::make_shared<W100FitFluxHandler>();
		}
		// W (111) case
		else if (materialName == "W111") {
			theFluxHandler = std::make_shared<W111FitFluxHandler>();
		}
		// Fe case
		else if (materialName == "Fe") {
			theFluxHandler = std::make_shared<FeFitFluxHandler>();
		}
		else {
			// The name was not recognized
			throw std::string(
					"\nThe name for the material is not known: \""
							+ materialName + "\"");
		}
	}

	else {
		std::cerr << "\nWarning: Material information has not been given."
				"\nDefaulting to tungsten with (100) surface \n" << std::endl;
		// we are to use a tungsten, W, flux handler
		theFluxHandler = std::make_shared<W100FitFluxHandler>();
	}

	// If the Helium fluence option is present, set the value
	if (options.useMaxHeliumFluence()) {
		theFluxHandler->setMaxHeFluence(options.getMaxHeliumFluence());
	}

	// If the helium flux option is present, set the value
	if (options.useHeliumFlux()) {
		theFluxHandler->setHeFlux(options.getHeliumFlux());
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

