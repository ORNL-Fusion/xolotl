#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class FeCrNetworkHandler : public NetworkHandler
{
public:
	FeCrNetworkHandler(const options::IOptions& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
