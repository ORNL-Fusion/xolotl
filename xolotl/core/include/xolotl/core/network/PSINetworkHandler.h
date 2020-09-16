#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class PSINetworkHandler : public NetworkHandler
{
public:
	PSINetworkHandler(const options::Options& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
