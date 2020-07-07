#pragma once

#include <xolotl/core/network/INetworkHandler.h>
#include <xolotl/options/Options.h>

namespace xolotl
{
namespace core
{
namespace network
{
using NetworkGeneratorFunction =
	std::function<std::shared_ptr<IReactionNetwork>(const options::Options&)>;

class NetworkHandler : public INetworkHandler
{
public:
	std::shared_ptr<IReactionNetwork>
	getNetwork() const final
	{
		return _network;
	}

protected:
	NetworkHandler(
		const options::Options& options, NetworkGeneratorFunction func);

protected:
	std::shared_ptr<IReactionNetwork> _network;
};
} // namespace network
} // namespace core
} // namespace xolotl
