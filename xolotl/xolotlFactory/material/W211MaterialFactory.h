#ifndef W211MATERIALHANDLERFACTORY_H
#define W211MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W211FitFluxHandler.h>
#include <W211AdvectionHandler.h>
#include <W211TrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (211) oriented tungsten material.
 */
class W211MaterialFactory: public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W211MaterialFactory(const xolotlCore::Options &options) :
		MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::W211FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W211AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W211TrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

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
