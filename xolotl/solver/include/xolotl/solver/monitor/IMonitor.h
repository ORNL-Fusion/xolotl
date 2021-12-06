#pragma once

namespace xolotl
{
namespace solver
{
namespace monitor
{
class IMonitor
{
public:
	virtual ~IMonitor()
	{
	}

	virtual void
	setup() = 0;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
