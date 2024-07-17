#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class T91NetworkHandler : public NetworkHandler
{
public:
	T91NetworkHandler(const options::IOptions& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
