#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace factory
{
namespace material
{
MaterialHandlerFactory::MaterialHandlerFactory() = default;

MaterialHandlerFactory::~MaterialHandlerFactory() = default;

MaterialHandlerFactory&
MaterialHandlerFactory::get()
{
	static MaterialHandlerFactory factory;
	return factory;
}

std::shared_ptr<core::material::IMaterialHandler>
MaterialHandlerFactory::generateMaterialHandler(const options::Options& options)
{
	const auto& materialName = options.getMaterial();
	auto it = _generators.find(materialName);
	if (it == _generators.end()) {
		throw std::runtime_error(
			"No material handler found for material \"" + materialName + "\"");
	}
	return it->second(options);
}

bool
MaterialHandlerFactory::registerGenerator(
	const std::string& name, const Generator& generator)
{
	return _generators.emplace(name, generator).second;
}
} // namespace material
} // namespace factory
} // namespace xolotl
