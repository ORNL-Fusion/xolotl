#ifndef OSTIMER_H
#define OSTIMER_H

#include <chrono>
#include <limits>

#include <xolotl/perf/ITimer.h>

namespace xolotl
{
namespace perf
{
namespace os
{
/**
 * A timer that measures how long something takes to execute.
 * Uses an operating system/runtime timer interface.
 */
class OSTimer : public ITimer
{
private:
	//! Concise name for type of our time source.
	using Clock = std::chrono::system_clock;

	//! Concise name for type of a timestamp.
	using Timestamp = std::chrono::time_point<Clock>;

	//! Concise name for type of a difference between timestamps.
	using Duration = std::chrono::duration<ITimer::ValType>;

	//! The timestamp value we use to indicate startTime is invalid.
	static const Timestamp invalidTimestamp;

	//! The timer's value.
	ITimer::ValType val;

	//! When the timer was started.
	//! Will be invalidTimestamp if timer is not running.
	Timestamp startTime;

public:
	OSTimer() : val(0), startTime(invalidTimestamp)
	{
	}

	virtual ~OSTimer(void)
	{
	}

	/**
	 * \see ITimer.h
	 */
	void
	start(void) override;

	/**
	 * \see ITimer.h
	 */
	void
	stop(void) override;

	/**
	 * \see ITimer.h
	 */
	void
	reset(void) override;

	/**
	 * Determine if the Timer is currently running.
	 *
	 * @return true if the Timer is running, false otherwise.
	 */
	virtual bool
	isRunning(void) const
	{
		return (startTime != invalidTimestamp);
	}

	/**
	 * \see ITimer.h
	 */
	ITimer::ValType
	getValue(void) const override
	{
		return val;
	}

	/**
	 * \see ITimer.h
	 */
	std::string
	getUnits(void) const override;

	/**
	 * Add the given Timer's value to my value.
	 *
	 * @param t The timer whose value should be added to my value.
	 * @return Myself after adding the given timer's value.
	 */
	virtual ITimer&
	operator+=(const ITimer& t)
	{
		val += t.getValue();
		return *this;
	}
};

} // namespace os
} // namespace perf
} // namespace xolotl

#endif // OSTIMER_H
