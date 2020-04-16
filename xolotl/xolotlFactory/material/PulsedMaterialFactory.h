#ifndef PULSEDMATERIALHANDLERFACTORY_H
#define PULSEDMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <PulsedFitFluxHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a tungsten material with
 * a intermittent incoming flux.
 */
class PulsedMaterialFactory: public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	PulsedMaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	PulsedMaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::PulsedFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

		// Switch on the dimension for the diffusion handler
		switch (dim) {
		case 0:
			theDiffusionHandler = std::make_shared<
					xolotlCore::DummyDiffusionHandler>();
			break;
		case 1:
			theDiffusionHandler = std::make_shared<
					xolotlCore::Diffusion1DHandler>();
			break;
		case 2:
			theDiffusionHandler = std::make_shared<
					xolotlCore::Diffusion2DHandler>();
			break;
		case 3:
			theDiffusionHandler = std::make_shared<
					xolotlCore::Diffusion3DHandler>();
			break;
		default:
			// The asked dimension is not good (e.g. -1, 4)
			throw std::string(
					"\nxolotlFactory: Bad dimension for the pulsed material factory.");
		}

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
