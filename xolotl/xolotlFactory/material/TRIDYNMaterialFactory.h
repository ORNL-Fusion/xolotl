#ifndef TRIDYNMATERIALHANDLERFACTORY_H
#define TRIDYNMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <TRIDYNFitFluxHandler.h>
#include <W100AdvectionHandler.h>
#include <W100TrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (100) oriented tungsten material with
 * a TRIDYN input file.
 */
class TRIDYNMaterialFactory: public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	TRIDYNMaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	TRIDYNMaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::TRIDYNFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W100AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W100TrapMutationHandler>();
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
					"\nxolotlFactory: Bad dimension for the TRIDYN material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~TRIDYNMaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // TRIDYNMATERIALHANDLERFACTORY_H
