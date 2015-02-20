#ifndef W111MATERIALHANDLERFACTORY_H
#define W111MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W111FitFluxHandler.h>
#include <W111AdvectionHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (111) oriented tungsten material.
 */
class W111MaterialFactory : public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W111MaterialFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W111MaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W111FitFluxHandler>();
		theAdvectionHandler = std::make_shared<xolotlCore::W111AdvectionHandler>();

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
				throw std::string("\nxolotlFactory: Bad dimension for the W111 material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W111MaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // W111MATERIALHANDLERFACTORY_H
