#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/viz/VizHandlerRegistryFactory.h>

namespace xolotl
{
namespace factory
{
using viz::VizHandlerRegistryFactory;
using xolotl::viz::IVizHandlerRegistry;

template VizHandlerRegistryFactory&
Factory<VizHandlerRegistryFactory, IVizHandlerRegistry>::get();

template std::shared_ptr<IVizHandlerRegistry>
Factory<VizHandlerRegistryFactory, IVizHandlerRegistry>::generate(
	const options::Options&);

template bool
Factory<VizHandlerRegistryFactory, IVizHandlerRegistry>::registerGenerator(
	const std::string&, const Generator&);
} // end namespace factory
} // end namespace xolotl
