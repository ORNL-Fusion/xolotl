#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class INetworkHandler;
}
} // namespace core

namespace factory
{
namespace network
{
class NetworkHandlerFactory :
	public Factory<NetworkHandlerFactory, core::network::INetworkHandler>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "NetworkHandlerFactory";
	}

	static std::string
	getName(const options::Options& options)
	{
		return options.getMaterial();
	}
};
} // namespace network
} // namespace factory
} // namespace xolotl
