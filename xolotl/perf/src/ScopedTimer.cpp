#include <xolotl/perf/ITimer.h>
#include <xolotl/perf/ScopedTimer.h>

namespace xolotl
{
namespace perf
{
ScopedTimer::ScopedTimer(std::shared_ptr<ITimer> _timer) : timer(_timer)
{
	timer->start();
}

ScopedTimer::~ScopedTimer()
{
	timer->stop();
}
} // end namespace perf
} // end namespace xolotl
