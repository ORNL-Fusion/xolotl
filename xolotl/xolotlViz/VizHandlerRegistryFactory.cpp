#include "XolotlConfigViz.h"
#include "VizHandlerRegistryFactory.h"
#include <DummyHandlerRegistry.h>
#include <iostream>

#if defined(HAVE_VIZLIB_STD)
#include <StandardHandlerRegistry.h>
#endif // defined(HAVE_VIZLIB_STD)

namespace xolotlViz {

static std::shared_ptr<IVizHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
bool initialize(bool useStdRegistry) {
	bool ret = true;

	if (useStdRegistry) {
#if defined(HAVE_VIZLIB_STD)
		// we are to use a standard handler registry
		theHandlerRegistry = std::make_shared<StandardHandlerRegistry>();
#else
		// it is not possible to use the standard registry
		throw std::string("\nxolotlViz::initialize: unable to build requested visualization standard handler registry due to missing dependencies");
#endif // defined(HAVE_VIZLIB_STD)
	} else {
		// we are to use a dummy handler registry
		theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
	}

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<IVizHandlerRegistry> getVizHandlerRegistry(void) {
	if (!theHandlerRegistry) {
		// Throw an error since we have not yet been initialized
		throw std::string("\nxolotlViz handler registry requested, but library has not been initialized");
	}

	return theHandlerRegistry;
}

} // end namespace xolotlViz
