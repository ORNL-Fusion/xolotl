#include <iostream>
#include <stdexcept>

#include <xolotl/perf/papi/PAPITimer.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
PAPITimer::Timestamp PAPITimer::invalidValue = -1;

PAPITimer::~PAPITimer(void)
{
	// nothing to do
}

void
PAPITimer::start(void)
{
	if (isRunning()) {
		throw std::runtime_error(
			"Attempting to start a timer that is already running.");
	}

	// Start the timer by sampling the current time.
	startTime = GetCurrentTime();
}

void
PAPITimer::stop(void)
{
	if (!isRunning()) {
		throw std::runtime_error(
			"Attempting to stop a timer that was not running.");
	}

	// Form the difference between the end timestamp and
	// our saved start timestamp.
	val += ToSeconds(GetCurrentTime() - startTime);

	// Indicate the timer is no longer running.
	startTime = invalidValue;
}

void
PAPITimer::reset(void)
{
	if (isRunning()) {
		throw std::runtime_error("Attempting to reset a timer that is running");
	}

	val = 0;
}

std::string
PAPITimer::getUnits(void) const
{
	return std::string("s");
}

} // namespace papi
} // namespace perf
} // namespace xolotl
