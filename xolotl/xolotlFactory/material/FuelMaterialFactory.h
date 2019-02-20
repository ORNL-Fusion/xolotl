#ifndef FUELMATERIALHANDLERFACTORY_H
#define FUELMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <FuelFitFluxHandler.h>
#include <ReSolutionHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a nuclear fuel material.
 */
class FuelMaterialFactory: public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	FuelMaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	FuelMaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::FuelFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theReSolutionHandler =
				std::make_shared<xolotlCore::ReSolutionHandler>();
		theDesorptionHandler = std::make_shared<
				xolotlCore::DummyDesorptionHandler>();

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
					"\nxolotlFactory: Bad dimension for the Fuel material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~FuelMaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 * \see IMaterialFactory.h
	 */
	void initializeMaterial(const xolotlCore::Options &options) {
		// Call the mother method
		MaterialFactory::initializeMaterial(options);

		// Change the flux amplitude because we have to take into account
		// that there are one xenon created every 4 fissions.
		theFluxHandler->setFluxAmplitude(options.getFluxAmplitude() / 4.0);

		return;
	}
};

} // end namespace xolotlFactory

#endif // FUELMATERIALHANDLERFACTORY_H
