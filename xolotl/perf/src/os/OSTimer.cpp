#include <stdexcept>

#include <xolotl/perf/os/OSTimer.h>

namespace xolotl
{
namespace perf
{
namespace os
{
const OSTimer::Timestamp OSTimer::invalidTimestamp =
	std::numeric_limits<Timestamp>::max();

void
OSTimer::start(void)
{
	if (isRunning()) {
		throw std::runtime_error(
			"Attempting to start a timer that is already running.");
	}

	// Start the timer by sampling the current time.
	startTime = Clock::now();
}

void
OSTimer::stop(void)
{
	if (!isRunning()) {
		throw std::runtime_error(
			"Attempting to stop a timer that was not running.");
	}

	// Stop the timer by sampling the ending time.
	auto endTime = Clock::now();

	// Form the difference between the end timestamp and
	// our saved start timestamp.
	val += static_cast<Duration>(endTime - startTime).count();

	// Indicate the timer is no longer running.
	startTime = invalidTimestamp;
}

void
OSTimer::reset(void)
{
	if (isRunning()) {
		throw std::runtime_error("Attempting to reset a timer that is running");
	}

	val = 0;
}

std::string
OSTimer::getUnits(void) const
{
	return std::string("s");
}

} // namespace os
} // namespace perf
} // namespace xolotl
