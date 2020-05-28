#ifndef PULSEDMATERIALHANDLERFACTORY_H
#define PULSEDMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/PulsedFitFluxHandler.h>
#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/modifiedreaction/trapmutation/DummyTrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a tungsten material with
 * a intermittent incoming flux.
 */
class PulsedMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	PulsedMaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::PulsedFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~PulsedMaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	void initializeMaterial(const xolotlCore::Options &options) {
		// Call the general method first
		MaterialFactory::initializeMaterial(options);

		// Set the pulse parameters
		theFluxHandler->setPulseTime(options.getPulseTime());
		theFluxHandler->setProportion(options.getPulseProportion());
	}
};

} // end namespace xolotlFactory

#endif // PULSEDMATERIALHANDLERFACTORY_H
