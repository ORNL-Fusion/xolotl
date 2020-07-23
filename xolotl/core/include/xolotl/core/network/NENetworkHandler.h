#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class NENetworkHandler : public NetworkHandler
{
public:
	NENetworkHandler(const options::Options& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
