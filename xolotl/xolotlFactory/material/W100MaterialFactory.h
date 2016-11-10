#ifndef W100MATERIALHANDLERFACTORY_H
#define W100MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W100FitFluxHandler.h>
#include <W100AdvectionHandler.h>
#include <W100TrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (100) oriented tungsten material.
 */
class W100MaterialFactory : public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W100MaterialFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W100MaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W100FitFluxHandler>();
		theAdvectionHandler.push_back(std::make_shared<xolotlCore::W100AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<xolotlCore::W100TrapMutationHandler>();

		// Switch on the dimension for the diffusion handler
		switch (dim) {
			case 1:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion1DHandler>();
				break;
			case 2:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion2DHandler>();
				break;
			case 3:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion3DHandler>();
				break;
			default:
				// The asked dimension is not good (e.g. -1, 0, 4)
				throw std::string("\nxolotlFactory: Bad dimension for the W100 material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W100MaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // W100MATERIALHANDLERFACTORY_H
