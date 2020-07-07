#include <iostream>

#include <xolotl/factory/viz/VizHandlerRegistryFactory.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace factory
{
namespace viz
{
VizHandlerRegistryFactory::VizHandlerRegistryFactory() = default;

VizHandlerRegistryFactory::~VizHandlerRegistryFactory() = default;

VizHandlerRegistryFactory&
VizHandlerRegistryFactory::get()
{
	static VizHandlerRegistryFactory factory;
	return factory;
}

std::shared_ptr<xolotl::viz::IVizHandlerRegistry>
VizHandlerRegistryFactory::generateVizHandlerRegistry(
	const options::Options& options)
{
	const std::string name = "";
	auto it = _generators.find(name);
	if (it == _generators.end()) {
		throw std::runtime_error("No viz handler found for \"" + name + "\"");
	}
	return it->second(options);
}

bool
VizHandlerRegistryFactory::registerGenerator(
	const std::string& name, const Generator& generator)
{
	return _generators.emplace(name, generator).second;
}
} // end namespace viz
} // end namespace factory
} // end namespace xolotl
