#ifndef FUELMATERIALHANDLERFACTORY_H
#define FUELMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <FuelFitFluxHandler.h>
#include <HeterogeneousNucleationHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a nuclear fuel material.
 */
class FuelMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	FuelMaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::FuelFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::HeterogeneousNucleationHandler>();

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
		theFluxHandler->setFluxAmplitude(
				options.getFluxAmplitude() * options.getFissionYield());

		// Pass the fission yield to the heterogenenous nucletation handler
		theNucleationHandler->setFissionYield(options.getFissionYield());

		return;
	}
};

} // end namespace xolotlFactory

#endif // FUELMATERIALHANDLERFACTORY_H
