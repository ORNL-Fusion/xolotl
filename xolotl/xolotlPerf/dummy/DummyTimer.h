#ifndef DUMMYTIMER_H
#define DUMMYTIMER_H

// Include
#include <string>
#include <memory>
#include "Timer.h"

using namespace std;

namespace xolotlPerf{

/**
 * The DummyTimer class is instantiated by the DummerHandlerRegistry class
 * and realizes the DummyTimer interface.
 */
class DummyTimer : public Timer
{
private:

	/** The defaule constructor is declared as private since Timers
	 *  must be initialized with a name.
	 */
	DummyTimer():Timer("") {}

public:

	/**
	 * DummyTimer constructor that takes the argument name
	 * to distinguish specific DummyTimer.
	 *
	 * @param aname The DummyTimer's name
	 */
	DummyTimer(std::string aname);
//	DummyTimer(std::string aname) : name(aname) {}

	/** The constructor.
	 */
	~DummyTimer();

	/**
	 * This operation returns a DummyTimer that is created using the copy
	 * constructor. If this DummyTimer is actually a subclass of DummyTimer, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * DummyTimer.
	 * @return A copy of this DummyTimer.
	 */
//	virtual std::shared_ptr<Timer> clone();

    // This operations starts the DummyTimer.
    void start();

    // This operation stops the DummyTimer.
    void stop();

	/**
	 * This operation returns the name of the Timer.
	 *
	 * @return The name of this Timer
	 */
	const std::string getName() const;

    // This operation returns the value of the Timer.
    double getValue();

    // This operation returns the units of the DummyTimer.
    long getUnits() const;

};  //end class DummyTimer

}  //end namespace xolotlPerf

#endif
