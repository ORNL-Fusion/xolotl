#ifndef ALLOYMATERIALHANDLERFACTORY_H
#define ALLOYMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/AlloyFitFluxHandler.h>
#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

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
	AlloyMaterialFactory(const options::Options &opts) :
			MaterialFactory(opts) {
		theFluxHandler = std::make_shared<core::flux::AlloyFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<core::advection::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				core::modified::DummyTrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~AlloyMaterialFactory() {
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // ALLOYMATERIALHANDLERFACTORY_H
