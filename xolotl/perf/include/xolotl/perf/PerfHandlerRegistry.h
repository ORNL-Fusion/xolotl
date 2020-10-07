#pragma once

#include <memory>

#include <xolotl/perf/IPerfHandler.h>

namespace xolotl
{
namespace perf
{
class PerfHandlerRegistry
{
public:
	static void
	set(const std::shared_ptr<IPerfHandler>& handler);

	static std::shared_ptr<IPerfHandler>
	get();

private:
	static std::shared_ptr<IPerfHandler> perfHandler;
};
} // namespace perf
} // namespace xolotl
