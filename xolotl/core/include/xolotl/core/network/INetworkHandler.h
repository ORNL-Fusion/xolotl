#pragma once

#include <memory>

#include <xolotl/core/network/IReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
class INetworkHandler
{
public:
	virtual ~INetworkHandler()
	{
	}

	virtual std::shared_ptr<IReactionNetwork>
	getNetwork() const = 0;
};

void
loadNetworkHandlers();
} // namespace network
} // namespace core
} // namespace xolotl
