#pragma once

#include <memory>

#include <xolotl/perf/ITimer.h>
#include <xolotl/util/Log.h>

namespace xolotl
{
namespace perf
{
class ITimer;

/**
 * A class for managing timer start/stop lifetime by code scope.
 * Used to simplify a common use case for a timer (starting timer when
 * enter a scope, and stopping timer when leave the scope regardless
 * of how we leave the scope).
 */
struct ScopedTimer
{
	/// The timer that should be active in the struct's scope.
	std::shared_ptr<ITimer> timer;

	ScopedTimer(std::shared_ptr<ITimer> _timer) : timer(_timer)
	{
		if (timer) {
			timer->start();
		}
		else {
			XOLOTL_LOG_XTRA << "Scoping a timer that does not exist";
		}
	}

	~ScopedTimer()
	{
		if (timer) {
			timer->stop();
		}
	}
};
} // end namespace perf
} // end namespace xolotl
