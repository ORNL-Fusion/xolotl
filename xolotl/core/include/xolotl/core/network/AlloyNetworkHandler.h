#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class AlloyNetworkHandler : public NetworkHandler
{
public:
	AlloyNetworkHandler(const options::Options& options);
};

namespace detail
{
auto alloyNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		AlloyNetworkHandler>({"800H"});
}
} // namespace network
} // namespace core
} // namespace xolotl
