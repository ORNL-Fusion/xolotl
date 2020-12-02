#pragma once

#include <xolotl/perf/ITimer.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
/**
 * A timer that measures how long something takes to execute.
 * Uses PAPI for sampling the system's timer.
 */
class PAPITimer : public ITimer
{
private:
	//! The type PAPI uses for timestamps.
	typedef long long Timestamp;

	//! An invalid Timestamp value.
	static constexpr Timestamp invalidValue = -1;

	//! The timer's value.
	ITimer::ValType val;

	//! When the timer was started.
	//! Will be invalidValue if timer is not running.
	Timestamp startTime;

	/**
	 * Sample the current time.
	 */
	Timestamp
	getCurrentTime() const;

	/**
	 * Convert a Timestamp to seconds.
	 */
	static ITimer::ValType
	toSeconds(Timestamp t)
	{
		// The timer we use - PAPI_get_real_nsec - gives time in nanoseconds.
		return t / 1.0e9;
	}

public:
	PAPITimer() : val(0), startTime(invalidValue)
	{
	}

	virtual ~PAPITimer()
	{
	}

	/**
	 * \see ITimer.h
	 */
	void
	start() override;

	/**
	 * \see ITimer.h
	 */
	void
	stop() override;

	/**
	 * \see ITimer.h
	 */
	void
	reset() override;

	/**
	 * Determine if the Timer is currently running.
	 *
	 * @return true if the Timer is running, false otherwise.
	 */
	bool
	isRunning() const
	{
		return (startTime != invalidValue);
	}

	/**
	 * \see ITimer.h
	 */
	ValType
	getValue() const override
	{
		return val;
	}

	/**
	 * \see ITimer.h
	 */
	std::string
	getUnits() const override;

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
} // namespace papi
} // namespace perf
} // namespace xolotl
