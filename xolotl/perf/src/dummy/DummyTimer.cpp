#include <xolotl/perf/dummy/DummyTimer.h>

namespace xolotl
{
namespace perf
{
namespace dummy
{
void
DummyTimer::start(void)
{
	// nothing to do in this stub routine
}

void
DummyTimer::stop(void)
{
	// nothing to do in this stub routine
}

void
DummyTimer::reset(void)
{
	// nothing to do in this stub routine
}

double
DummyTimer::getValue(void) const
{
	return 0;
}

std::string
DummyTimer::getUnits(void) const
{
	return "";
}

} // namespace dummy
} // namespace perf
} // namespace xolotl
