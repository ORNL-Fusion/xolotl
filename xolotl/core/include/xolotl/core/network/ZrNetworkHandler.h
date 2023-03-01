#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class ZrNetworkHandler : public NetworkHandler
{
public:
	ZrNetworkHandler(const options::IOptions& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
