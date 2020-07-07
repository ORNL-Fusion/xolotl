#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace factory
{
namespace network
{
NetworkHandlerFactory::NetworkHandlerFactory() = default;

NetworkHandlerFactory::~NetworkHandlerFactory() = default;

NetworkHandlerFactory&
NetworkHandlerFactory::get()
{
	static NetworkHandlerFactory factory;
	return factory;
}

std::shared_ptr<core::network::INetworkHandler>
NetworkHandlerFactory::generateNetworkHandler(const options::Options& options)
{
	const auto& materialName = options.getMaterial();
	auto it = _generators.find(materialName);
	if (it == _generators.end()) {
		throw std::runtime_error(
			"No network handler found for material \"" + materialName + "\"");
	}
	return it->second(options);
}

bool
NetworkHandlerFactory::registerGenerator(
	const std::string& name, const Generator& generator)
{
	return _generators.emplace(name, generator).second;
}
} // namespace network
} // namespace factory
} // namespace xolotl
