#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class LiNetworkHandler : public NetworkHandler
{
public:
	LiNetworkHandler(const options::IOptions& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
