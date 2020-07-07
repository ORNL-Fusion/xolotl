#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class FeNetworkHandler : public NetworkHandler
{
public:
	FeNetworkHandler(const options::Options& options);
};

namespace detail
{
auto feNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		FeNetworkHandler>({"Fe"});
}
} // namespace network
} // namespace core
} // namespace xolotl
