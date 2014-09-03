#include "MaterialHandlerFactory.h"
#include "W100FitFluxHandler.h"
#include "W110FitFluxHandler.h"
#include "W111FitFluxHandler.h"
#include "W210FitFluxHandler.h"
#include "W211FitFluxHandler.h"
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
		// W (110) case
		else if (materialName == "W110") {
			std::cout << "\nThe material being used is W with a (110) surface orientation." << std::endl;
			theFluxHandler = std::make_shared<W110FitFluxHandler>();
		}
		// W (111) case
		else if (materialName == "W111") {
			theFluxHandler = std::make_shared<W111FitFluxHandler>();
		}
		// W (210) case
		else if (materialName == "W210") {
			std::cout << "\nThe material: " << materialName << std::endl;
			theFluxHandler = std::make_shared<W210FitFluxHandler>();
		}
		// W (211) case
		else if (materialName == "W211") {
			std::cout << "\nThe material: " << materialName << std::endl;
			theFluxHandler = std::make_shared<W211FitFluxHandler>();
		}
		// Fe case
		else if (materialName == "Fe") {
			theFluxHandler = std::make_shared<FeFitFluxHandler>();
		}
		else {
			// The name was not recognized
			throw std::string(
					"\nThe name material name is not known: \""
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

