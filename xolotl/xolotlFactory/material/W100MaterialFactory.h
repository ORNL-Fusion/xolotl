#ifndef W100MATERIALHANDLERFACTORY_H
#define W100MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W100FitFluxHandler.h>
#include <W100AdvectionHandler.h>
#include <DiffusionHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (100) oriented tungsten material.
 */
class W100MaterialFactory : public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 */
	W100MaterialFactory() {
		theFluxHandler = std::make_shared<xolotlCore::W100FitFluxHandler>();
		theAdvectionHandler = std::make_shared<xolotlCore::W100AdvectionHandler>();
		theDiffusionHandler = std::make_shared<xolotlCore::DiffusionHandler>();
		return;
	}

	/**
	 * The destructor
	 */
	~W100MaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // W100MATERIALHANDLERFACTORY_H
