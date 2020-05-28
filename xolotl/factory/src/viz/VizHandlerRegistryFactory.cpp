#include <iostream>
#include <mpi.h>
#include <xolotl/viz/config.h>
#include <xolotl/factory/viz/VizHandlerRegistryFactory.h>
#include <xolotl/viz/dummy/DummyHandlerRegistry.h>

#if defined(HAVE_VIZLIB_STD)
#include <xolotl/viz/standard/StandardHandlerRegistry.h>
#endif // defined(HAVE_VIZLIB_STD)

namespace xolotlFactory {

static std::shared_ptr<xolotlViz::IVizHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
bool initializeVizHandler(bool useStdRegistry) {
	bool ret = true;

	if (useStdRegistry) {
#if defined(HAVE_VIZLIB_STD)
		// we are to use a standard handler registry
		theHandlerRegistry = std::make_shared<xolotlViz::StandardHandlerRegistry>();
#else
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		// Only print the error message once when running in parallel
		if (procId == 0) {
			// it is not possible to use the standard registry
			throw std::string(
					"\nxolotlFactory::initialize: unable to build requested visualization "
							"standard handler registry due to missing dependencies");
		}
#endif // defined(HAVE_VIZLIB_STD)
	} else {
		// we are to use a dummy handler registry
		theHandlerRegistry = std::make_shared<xolotlViz::DummyHandlerRegistry>();
	}

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<xolotlViz::IVizHandlerRegistry> getVizHandlerRegistry(void) {
	if (!theHandlerRegistry) {
		// Throw an error since we have not yet been initialized
		throw std::string(
				"\nxolotlFactory handler registry requested, but library has not been initialized.");
	}

	return theHandlerRegistry;
}

} // end namespace xolotlFactory
