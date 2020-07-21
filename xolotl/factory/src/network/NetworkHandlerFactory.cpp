#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace factory
{
using core::network::INetworkHandler;
using network::NetworkHandlerFactory;

template NetworkHandlerFactory&
Factory<NetworkHandlerFactory, INetworkHandler>::get();

template std::shared_ptr<INetworkHandler>
Factory<NetworkHandlerFactory, INetworkHandler>::generate(
	const options::Options&);

template bool
Factory<NetworkHandlerFactory, INetworkHandler>::registerGenerator(
	const std::string&, const Generator&);
} // namespace factory
} // namespace xolotl
