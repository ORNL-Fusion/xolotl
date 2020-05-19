#ifndef ALLOYMATERIALHANDLERFACTORY_H
#define ALLOYMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <AlloyFitFluxHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for an alloy material.
 */
class AlloyMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	AlloyMaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::AlloyFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

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
