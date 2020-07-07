#include <xolotl/util/MPIUtils.h>
#include <xolotl/viz/VizHandlerRegistry.h>
#include <xolotl/viz/config.h>
#include <xolotl/viz/dummy/DummyHandlerRegistry.h>

#if defined(HAVE_VIZLIB_STD)
#include <xolotl/viz/standard/StandardHandlerRegistry.h>
#endif // defined(HAVE_VIZLIB_STD)

namespace xolotl
{
namespace viz
{
std::shared_ptr<IVizHandlerRegistry>
	VizHandlerRegistry::staticVizHandlerRegistry;

VizHandlerRegistry::VizHandlerRegistry(const options::Options& options)
{
	if (options.useVizStandardHandlers()) {
#if defined(HAVE_VIZLIB_STD)
		// we are to use a standard handler registry
		staticVizHandlerRegistry =
			std::make_shared<standard::StandardHandlerRegistry>();
#else
		// Get the current process ID
		int procId;
		auto xolotlComm = util::getMPIComm();
		MPI_Comm_rank(xolotlComm, &procId);
		// Only print the error message once when running in parallel
		if (procId == 0) {
			// it is not possible to use the standard registry
			throw std::string(
				"\nxolotlFactory::initialize: unable to build requested "
				"visualization "
				"standard handler registry due to missing dependencies");
		}
#endif // defined(HAVE_VIZLIB_STD)
	}
	else {
		// we are to use a dummy handler registry
		staticVizHandlerRegistry =
			std::make_shared<dummy::DummyHandlerRegistry>();
	}
}

std::shared_ptr<IPlot>
VizHandlerRegistry::getPlot(const std::string& name, PlotType type)
{
	return staticVizHandlerRegistry->getPlot(name, type);
}

namespace detail
{
auto defaultVizHandlerRegistryRegistration =
	xolotl::factory::viz::VizHandlerRegistryFactory::Registration<
		VizHandlerRegistry>("");
}
} // namespace viz
} // namespace xolotl
