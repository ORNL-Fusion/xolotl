#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace perf
{
class IPerfHandler;
}

namespace factory
{
namespace perf
{
class PerfHandlerFactory :
	public Factory<PerfHandlerFactory, xolotl::perf::IPerfHandler>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "PerfHandlerFactory";
	}

	static std::string
	getName(const options::IOptions& options)
	{
		return options.getPerfHandlerName();
	}
};
} // namespace perf
} // namespace factory
} // namespace xolotl
