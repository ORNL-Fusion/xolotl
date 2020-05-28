#ifndef DUMMYTIMER_H
#define DUMMYTIMER_H

#include <xolotl/perf/ITimer.h>
#include <xolotl/core/Identifiable.h>

using namespace std;

namespace xolotlPerf {

/**
 * The DummyTimer class is instantiated by the DummerHandlerRegistry class
 * and realizes the DummyTimer interface.
 */
class DummyTimer: public ITimer, public xolotlCore::Identifiable {
private:

	/**
	 * The default constructor is declared as private since Timers
	 *  must be initialized with a name.
	 */
	DummyTimer(void) :
			xolotlCore::Identifiable("unused") {
	}

public:

	/**
	 * DummyTimer constructor that takes the argument timerName
	 * to distinguish specific DummyTimer.
	 *
	 * @param name The DummyTimer's name
	 */
	DummyTimer(const std::string& name) :
			xolotlCore::Identifiable("unused") {
	}

	/**
	 * Destroy the timer.
	 */
	virtual ~DummyTimer(void) {
	}

	/**
	 * Start the timer.
	 */
	virtual void start(void);

	/**
	 * Stop the timer.
	 */
	virtual void stop(void);

	/**
	 * Reset the timer's value.
	 */
	virtual void reset(void);

	/**
	 * Obtain the timer's value.
	 */
	virtual ITimer::ValType getValue(void) const;

	/**
	 * Obtain a string describing the units of the timer's value.
	 */
	virtual std::string getUnits(void) const;

};
//end class DummyTimer

}//end namespace xolotlPerf

#endif
