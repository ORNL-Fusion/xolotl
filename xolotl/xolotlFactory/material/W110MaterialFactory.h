#ifndef W110MATERIALHANDLERFACTORY_H
#define W110MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W110FitFluxHandler.h>
#include <W110AdvectionHandler.h>
#include <DiffusionHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (110) oriented tungsten material.
 */
class W110MaterialFactory : public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 */
	W110MaterialFactory() {
		theFluxHandler = std::make_shared<xolotlCore::W110FitFluxHandler>();
		theAdvectionHandler = std::make_shared<xolotlCore::W110AdvectionHandler>();
		theDiffusionHandler = std::make_shared<xolotlCore::DiffusionHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W110MaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // W110MATERIALHANDLERFACTORY_H
