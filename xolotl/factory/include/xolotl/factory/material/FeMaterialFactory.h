#ifndef FEMATERIALHANDLERFACTORY_H
#define FEMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/FeFitFluxHandler.h>
#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

/**
 * Subclass of MaterialFactory for an iron material.
 */
class FeMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	FeMaterialFactory(const options::Options &opts) :
			MaterialFactory(opts) {
		theFluxHandler = std::make_shared<core::flux::FeFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<core::advection::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				core::modified::DummyTrapMutationHandler>();
		theNucleationHandler = std::make_shared<
				core::modified::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~FeMaterialFactory() {
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // FEMATERIALHANDLERFACTORY_H
