#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class RPVNetworkHandler : public NetworkHandler
{
public:
	RPVNetworkHandler(const options::IOptions& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
