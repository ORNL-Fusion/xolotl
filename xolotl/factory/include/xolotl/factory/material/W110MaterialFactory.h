#ifndef W110MATERIALHANDLERFACTORY_H
#define W110MATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/W110FitFluxHandler.h>
#include <xolotl/core/advection/W110AdvectionHandler.h>
#include <xolotl/core/modifiedreaction/trapmutation/W110TrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (110) oriented tungsten material.
 */
class W110MaterialFactory: public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W110MaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::W110FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W110AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W110TrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W110MaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // W110MATERIALHANDLERFACTORY_H
