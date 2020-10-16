#pragma once

#include <xolotl/perf/PerfHandler.h>

namespace xolotl
{
namespace perf
{
namespace os
{
class OSHandler : public PerfHandler
{
public:
	using PerfHandler::PerfHandler;

	virtual ~OSHandler();

	std::shared_ptr<ITimer>
	getTimer(const std::string& name) override;
};
} // namespace os
} // namespace perf
} // namespace xolotl
