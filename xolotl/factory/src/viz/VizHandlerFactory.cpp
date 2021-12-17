#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/viz/VizHandlerFactory.h>

namespace xolotl
{
namespace factory
{
using viz::VizHandlerFactory;
using xolotl::viz::IVizHandler;

template VizHandlerFactory&
Factory<VizHandlerFactory, IVizHandler>::get(const UserInitializer&);

template std::shared_ptr<IVizHandler>
Factory<VizHandlerFactory, IVizHandler>::generate(const options::IOptions&);

template bool
Factory<VizHandlerFactory, IVizHandler>::registerGenerator(
	const std::string&, const Generator&);
} // end namespace factory
} // end namespace xolotl
