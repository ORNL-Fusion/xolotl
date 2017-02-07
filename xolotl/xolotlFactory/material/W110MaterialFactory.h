#ifndef W110MATERIALHANDLERFACTORY_H
#define W110MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W110FitFluxHandler.h>
#include <W110AdvectionHandler.h>
#include <W110TrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (110) oriented tungsten material.
 */
class W110MaterialFactory : public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W110MaterialFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W110MaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W110FitFluxHandler>();
		theAdvectionHandler.push_back(std::make_shared<xolotlCore::W110AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<xolotlCore::W110TrapMutationHandler>();

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
				throw std::string("\nxolotlFactory: Bad dimension for the W110 material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W110MaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // W110MATERIALHANDLERFACTORY_H
