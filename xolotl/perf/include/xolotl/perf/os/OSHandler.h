#pragma once

#include <xolotl/perf/standard/StdHandler.h>

namespace xolotl
{
namespace perf
{
namespace os
{
class OSHandler : public standard::StdHandler
{
public:
	OSHandler(const options::IOptions&);

	virtual ~OSHandler()
	{
	}

	std::shared_ptr<ITimer>
	getTimer(const std::string& name) override;
};
} // namespace os
} // namespace perf
} // namespace xolotl
