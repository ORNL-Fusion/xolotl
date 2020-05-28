#ifndef PAPITIMER_H
#define PAPITIMER_H

#include <xolotl/perf/config.h>
#if !defined(HAVE_PAPI)
#  error "Using PAPI-based handler registry classes but PAPI was not found when configured."
#endif // !defined(HAVE_PAPI)

#include <xolotl/perf/ITimer.h>
#include <xolotl/core/Identifiable.h>
#include <papi.h>

namespace xolotlPerf {

/// A timer that measures how long something takes to execute.
/// Uses PAPI for sampling the system's timer.
class PAPITimer: public ITimer, public xolotlCore::Identifiable {
private:
	/// The type PAPI uses for timestamps.
	typedef long long Timestamp;

	/// An invalid Timestamp value.
	static Timestamp invalidValue;

	/// The timer's value.
	ITimer::ValType val;

	/// When the timer was started.
	/// Will be invalidValue if timer is not running.
	Timestamp startTime;

	/// Construct a timer.
	/// The default constructor is private to force callers to provide a name for the timer object.
	PAPITimer(void) :
			xolotlCore::Identifiable("unused"), val(0) {
	}

	/// Sample the current time.
	///
	/// @return The current time.
	Timestamp GetCurrentTime(void) const {
		return PAPI_get_real_nsec();
	}

	/// Convert a Timestamp to seconds.
	///
	/// @param t A timestamp to be converted to seconds.
	/// @return t in terms of seconds.
	static ITimer::ValType ToSeconds(Timestamp t) {
		// The timer we use - PAPI_get_real_nsec - gives time in nanoseconds.
		return t / 1.0e9;
	}

public:
	///
	/// Construct a timer.
	///
	/// @param name The name to associate with the timer.
	PAPITimer(const std::string& name) :
			xolotlCore::Identifiable(name), val(0), startTime(invalidValue) {
	}

	///
	/// Destroy the timer.
	///
	virtual ~PAPITimer(void);

	///
	/// Start the timer.
	/// Throws std::runtime_error if starting a timer that was already started.
	///
	virtual void start(void);

	///
	/// Stop the timer.
	/// Throws std::runtime_error if stopping a timer that was not running.
	///
	virtual void stop(void);

	///
	/// Reset the timer.
	/// Throws std::runtime_error if timer is running.
	///
	virtual void reset(void);

	///
	/// Determine if the Timer is currently running.
	///
	/// @return true if the Timer is running, false otherwise.
	///
	virtual bool isRunning(void) const {
		return (startTime != invalidValue);
	}

	///
	/// Retrieve the value of the timer.
	/// The value is only valid if the timer is not running.
	///
	/// @return The elapsed time measured by this timer.
	///
	virtual ValType getValue(void) const {
		return val;
	}

	///
	/// Retrieve the Timer value's units.
	/// @return The units in which the timer's value is given.
	///
	virtual std::string getUnits(void) const;

	/// Add the given Timer's value to my value.
	/// @param t The timer whose value should be added to my value.
	/// @return Myself after adding the given timer's value.
	virtual ITimer& operator+=(const ITimer& t) {
		val += t.getValue();
		return *this;
	}
};

} // namespace xolotlPerf

#endif // PAPITIMER_H
