#ifndef IHANDLERREGISTRY_H
#define IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "Timer.h" //Dependency Generated Source:IHandlerRegistry Target:Timer
#include "EventCounter.h" //Dependency Generated Source:IHandlerRegistry Target:EventCounter
#include "HardwareCounter.h" //Dependency Generated Source:IHandlerRegistry Target:HardwareCounter
#include "HardwareQuantities.h"

namespace xolotlPerf {

class Timer;
class EventCounter;
class HardwareCounter;

// Realizations of this interface are responsible for the collection of performance data.
//class IHandlerRegistry //: public HardwareQuantities
//class IHandlerRegistry: public Timer, EventCounter, HardwareCounter {
class IHandlerRegistry {

public:

	/**
	 * The copy constructor
	 * @param other The HandlerRegistry to copy
	 */
	IHandlerRegistry(const IHandlerRegistry &other);

	// The destructor
	virtual ~IHandlerRegistry();

	// This operation returns the Timer specified by the parameter.
	virtual std::shared_ptr<Timer> getTimer(std::string name) const = 0;

	// This operation returns the EventCounter specified by the parameter.
	virtual std::shared_ptr<EventCounter> getEventCounter(
			std::string name) const = 0;

	// This operation returns the specified HardwareCounter.
	virtual std::shared_ptr<HardwareCounter> getHardwareCounter(
			std::string name,
			std::vector<HardwareQuantities> quantities) const = 0;

	// This operation returns a list of values of the, initially specified, PAPI preset quantities monitored by the HardwareCounter.
	virtual std::vector<HardwareQuantities> getHardwareQuantities() const = 0;

	// This operation outputs the information gathered.
	virtual void dump(std::ostream out) const = 0;

};
//end class IHandlerRegistry

}//end namespace xolotlPerff

#endif
