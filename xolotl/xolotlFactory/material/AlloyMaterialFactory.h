#ifndef ALLOYMATERIALHANDLERFACTORY_H
#define ALLOYMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <AlloyFitFluxHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for an alloy material.
 */
class AlloyMaterialFactory: public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	AlloyMaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	AlloyMaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::AlloyFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theReSolutionHandler = std::make_shared<
				xolotlCore::DummyReSolutionHandler>();
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
					"\nxolotlFactory: Bad dimension for the Alloy material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~AlloyMaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // ALLOYMATERIALHANDLERFACTORY_H
