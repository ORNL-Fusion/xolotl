#ifndef PULSEDMATERIALHANDLERFACTORY_H
#define PULSEDMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/MaterialFactory.h>
#include <xolotl/core/flux/PulsedFitFluxHandler.h>
#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

/**
 * Subclass of MaterialFactory for a tungsten material with
 * a intermittent incoming flux.
 */
class PulsedMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	PulsedMaterialFactory(const options::Options &opts) :
			MaterialFactory(opts) {
		theFluxHandler = std::make_shared<core::flux::PulsedFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<core::advection::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				core::modified::DummyTrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~PulsedMaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	void initializeMaterial(const options::Options &opts) {
		// Call the general method first
		MaterialFactory::initializeMaterial(opts);

		// Set the pulse parameters
		theFluxHandler->setPulseTime(opts.getPulseTime());
		theFluxHandler->setProportion(opts.getPulseProportion());
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // PULSEDMATERIALHANDLERFACTORY_H
