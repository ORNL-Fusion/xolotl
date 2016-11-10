#ifndef EVENTCOUNTER_H
#define EVENTCOUNTER_H

#include "xolotlCore/Identifiable.h"
#include "xolotlPerf/IEventCounter.h"

namespace xolotlPerf {

/**
 * An EventCounter keeps a count.  Code using an EventCounter can
 * increment the counter whenever an event of interest occurs, and
 * retrieve the count whenever it is desired.
 */
class EventCounter: public IEventCounter, public xolotlCore::Identifiable {
private:

	/**
	 * The value of this IEventCounter.
	 */
	IEventCounter::ValType value;

	/**
	 * We declare a private default constructor to force
	 * client code to provide a name when creating EventCounters.
	 */
	EventCounter(void) :
			xolotlCore::Identifiable("unused"), value(0) {
	}

public:

	/**
	 * EventCounter constructor that takes the argument name
	 *
	 * @param name The EventCounter's name
	 */
	EventCounter(const std::string& name) :
			xolotlCore::Identifiable(name), value(0) {
	}

	/**
	 * The destructor
	 */
	virtual ~EventCounter() {
	}

	/**
	 * This operation returns the value of the EventCounter,
	 * the frequency of the specified event.
	 */
	virtual IEventCounter::ValType getValue() const {
		return value;
	}

	/**
	 * This operation increments the EventCounter.
	 */
	virtual void increment() {
		++value;
	}

};
//end class EventCounter

}//end namespace xolotlPerf

#endif
