#include "MaterialHandlerFactory.h"
#include "W100FitFluxHandler.h"
#include "W110FitFluxHandler.h"
#include "W111FitFluxHandler.h"
#include "W210FitFluxHandler.h"
#include "W211FitFluxHandler.h"
#include "W221FitFluxHandler.h"
#include "W310FitFluxHandler.h"
#include "W311FitFluxHandler.h"
#include "W320FitFluxHandler.h"
#include "W321FitFluxHandler.h"
#include "FeFitFluxHandler.h"
#include <FluxHandler.h>
#include <iostream>

namespace xolotlSolver {

static std::shared_ptr<IFluxHandler> theFluxHandler;

bool initializeMaterial(xolotlCore::Options &options) {
	bool ret = true;

	if (options.useMaterial()) {
		std::string materialName = options.getMaterial();
		std::cout << "\nThe material being used is: " << materialName << std::endl;

		if ( materialName.substr(0,1) == "W" ) {

			std::string surfaceOrientation = materialName.substr(1, 3);
			int orientation = std::stoi(surfaceOrientation);

			switch ( orientation ) {
			case 100:
				theFluxHandler = std::make_shared<W100FitFluxHandler>();
				break;
			case 110:
				theFluxHandler = std::make_shared<W110FitFluxHandler>();
				break;
			case 111:
				theFluxHandler = std::make_shared<W111FitFluxHandler>();
				break;
			case 210:
				theFluxHandler = std::make_shared<W210FitFluxHandler>();
				break;
			case 211:
				theFluxHandler = std::make_shared<W211FitFluxHandler>();
				break;
			case 221:
				theFluxHandler = std::make_shared<W221FitFluxHandler>();
				break;
			case 310:
				theFluxHandler = std::make_shared<W310FitFluxHandler>();
				break;
			case 311:
				theFluxHandler = std::make_shared<W311FitFluxHandler>();
				break;
			case 320:
				theFluxHandler = std::make_shared<W320FitFluxHandler>();
				break;
			case 321:
				theFluxHandler = std::make_shared<W321FitFluxHandler>();
				break;

			}

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

