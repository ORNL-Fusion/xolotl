#ifndef W111MATERIALHANDLERFACTORY_H
#define W111MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W111FitFluxHandler.h>
#include <W111AdvectionHandler.h>
#include <W111TrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (111) oriented tungsten material.
 */
class W111MaterialFactory: public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W111MaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::W111FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W111AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W111TrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W111MaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // W111MATERIALHANDLERFACTORY_H
