#include <papi.h>

#include <iostream>
#include <stdexcept>

#include <xolotl/perf/papi/PAPITimer.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
PAPITimer::Timestamp
PAPITimer::getCurrentTime() const
{
	return PAPI_get_real_nsec();
}

void
PAPITimer::start()
{
	if (isRunning()) {
		throw std::runtime_error(
			"Attempting to start a timer that is already running.");
	}

	// Start the timer by sampling the current time.
	startTime = getCurrentTime();
}

void
PAPITimer::stop()
{
	if (!isRunning()) {
		throw std::runtime_error(
			"Attempting to stop a timer that was not running.");
	}

	// Form the difference between the end timestamp and
	// our saved start timestamp.
	val += toSeconds(getCurrentTime() - startTime);

	// Indicate the timer is no longer running.
	startTime = invalidValue;
}

void
PAPITimer::reset()
{
	if (isRunning()) {
		throw std::runtime_error("Attempting to reset a timer that is running");
	}

	val = 0;
}

std::string
PAPITimer::getUnits() const
{
	return std::string("s");
}

} // namespace papi
} // namespace perf
} // namespace xolotl
