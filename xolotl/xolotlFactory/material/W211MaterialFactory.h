#ifndef W211MATERIALHANDLERFACTORY_H
#define W211MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W211FitFluxHandler.h>
#include <W211AdvectionHandler.h>
#include <W211TrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (211) oriented tungsten material.
 */
class W211MaterialFactory: public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W211MaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W211MaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W211FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W211AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W211TrapMutationHandler>();
		theReSolutionHandler = std::make_shared<xolotlCore::DummyReSolutionHandler>();

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
					"\nxolotlFactory: Bad dimension for the W211 material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W211MaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // W211MATERIALHANDLERFACTORY_H
