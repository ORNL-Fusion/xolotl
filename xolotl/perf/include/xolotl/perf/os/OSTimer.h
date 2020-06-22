#ifndef OSTIMER_H
#define OSTIMER_H

#include <chrono>
#include <limits>

#include <xolotl/perf/ITimer.h>
#include <xolotl/perf/config.h>
#include <xolotl/util/Identifiable.h>

namespace xolotl
{
namespace perf
{
namespace os
{
/// A timer that measures how long something takes to execute.
/// Uses an operating system/runtime timer interface.
class OSTimer : public ITimer, public util::Identifiable
{
private:
	/// Concise name for type of our time source.
	using Clock = std::chrono::system_clock;

	/// Concise name for type of a timestamp.
	using Timestamp = std::chrono::time_point<Clock>;

	/// Concise name for type of a difference between timestamps.
	using Duration = std::chrono::duration<ITimer::ValType>;

	/// The timestamp value we use to indicate startTime is invalid.
	static const Timestamp invalidTimestamp;

	/// The timer's value.
	ITimer::ValType val;

	/// When the timer was started.
	/// Will be invalidTimestamp if timer is not running.
	Timestamp startTime;

	/// Construct a timer.
	/// The default constructor is private to force callers to provide a name
	/// for the timer object.
	OSTimer(void) : util::Identifiable("unused"), val(0)
	{
	}

public:
	///
	/// Construct a timer.
	///
	/// @param name The name to associate with the timer.
	OSTimer(const std::string& name) :
		util::Identifiable(name), val(0), startTime(invalidTimestamp)
	{
	}

	///
	/// Destroy the timer.
	///
	virtual ~OSTimer(void)
	{
	}

	///
	/// Start the timer.
	/// Throws std::runtime_error if starting a timer that was already started.
	///
	void
	start(void) override;

	///
	/// Stop the timer.
	/// Throws std::runtime_error if stopping a timer that was not running.
	///
	void
	stop(void) override;

	///
	/// Reset the timer's value.
	/// Throws std::runtime_error if resetting a timer that was running.
	///
	void
	reset(void) override;

	///
	/// Determine if the Timer is currently running.
	///
	/// @return true if the Timer is running, false otherwise.
	///
	virtual bool
	isRunning(void) const
	{
		return (startTime != invalidTimestamp);
	}
	///
	/// Retrieve the value of the timer.
	/// The value is only valid if the timer is not running.
	///
	/// @return The elapsed time measured by this timer.
	///
	ITimer::ValType
	getValue(void) const override
	{
		return val;
	}

	///
	/// Retrieve the Timer value's units.
	/// @return The units in which the timer's value is given.
	///
	std::string
	getUnits(void) const override;

	/// Add the given Timer's value to my value.
	/// @param t The timer whose value should be added to my value.
	/// @return Myself after adding the given timer's value.
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
