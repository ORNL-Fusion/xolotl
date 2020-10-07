#pragma once

#include <xolotl/perf/ITimer.h>
#include <xolotl/util/Identifiable.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
/// A timer that measures how long something takes to execute.
/// Uses PAPI for sampling the system's timer.
class PAPITimer : public ITimer, public util::Identifiable
{
private:
	/// The type PAPI uses for timestamps.
	typedef long long Timestamp;

	/// An invalid Timestamp value.
	static constexpr Timestamp invalidValue = -1;

	/// The timer's value.
	ITimer::ValType val;

	/// When the timer was started.
	/// Will be invalidValue if timer is not running.
	Timestamp startTime;

	/// Construct a timer.
	/// The default constructor is private to force callers to provide a name
	/// for the timer object.
	PAPITimer() : util::Identifiable("unused"), val(0)
	{
	}

	/// Sample the current time.
	///
	/// @return The current time.
	Timestamp
	getCurrentTime() const;

	/// Convert a Timestamp to seconds.
	///
	/// @param t A timestamp to be converted to seconds.
	/// @return t in terms of seconds.
	static ITimer::ValType
	toSeconds(Timestamp t)
	{
		// The timer we use - PAPI_get_real_nsec - gives time in nanoseconds.
		return t / 1.0e9;
	}

public:
	///
	/// Construct a timer.
	///
	/// @param name The name to associate with the timer.
	PAPITimer(const std::string& name);

	///
	/// Destroy the timer.
	///
	virtual ~PAPITimer();

	///
	/// Start the timer.
	/// Throws std::runtime_error if starting a timer that was already started.
	///
	void
	start() override;

	///
	/// Stop the timer.
	/// Throws std::runtime_error if stopping a timer that was not running.
	///
	void
	stop() override;

	///
	/// Reset the timer.
	/// Throws std::runtime_error if timer is running.
	///
	void
	reset() override;

	///
	/// Determine if the Timer is currently running.
	///
	/// @return true if the Timer is running, false otherwise.
	///
	bool
	isRunning() const
	{
		return (startTime != invalidValue);
	}

	///
	/// Retrieve the value of the timer.
	/// The value is only valid if the timer is not running.
	///
	/// @return The elapsed time measured by this timer.
	///
	ValType
	getValue() const override
	{
		return val;
	}

	///
	/// Retrieve the Timer value's units.
	/// @return The units in which the timer's value is given.
	///
	std::string
	getUnits() const override;

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
} // namespace papi
} // namespace perf
} // namespace xolotl
