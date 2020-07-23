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
} // namespace network
} // namespace core
} // namespace xolotl
