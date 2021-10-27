#pragma once

#include <memory>

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

	ScopedTimer(std::shared_ptr<ITimer> _timer);

	~ScopedTimer();
};
} // end namespace perf
} // end namespace xolotl
