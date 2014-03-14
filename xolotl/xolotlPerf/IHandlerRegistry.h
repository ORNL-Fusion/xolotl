#ifndef IHANDLERREGISTRY_H
#define IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "ITimer.h" //Dependency Generated Source:IHandlerRegistry Target:ITimer
#include "IEventCounter.h" //Dependency Generated Source:IHandlerRegistry Target:IEventCounter
#include "IHardwareCounter.h" //Dependency Generated Source:IHandlerRegistry Target:IHardwareCounter
#include "HardwareQuantities.h"


namespace xolotlPerf {

/**
 * Factory for building performance data collection objects, such
 * as timers and counters.
 */
class IHandlerRegistry {

public:

	/**
	 * The destructor
	 */
	virtual ~IHandlerRegistry(){}

	/**
	 * This operation returns the ITimer specified by the parameter.
	 */
	virtual std::shared_ptr<ITimer> getTimer(std::string name) = 0;

	/**
	 * This operation returns the IEventCounter specified by the parameter.
	 */
	virtual std::shared_ptr<IEventCounter> getEventCounter( std::string name) = 0;

	/**
	 * This operation returns the specified IHardwareCounter.
	 */
	virtual std::shared_ptr<IHardwareCounter> getHardwareCounter( std::string name,
			            std::vector<HardwareQuantities> quantities) = 0;

	/**
	 * This operation outputs the information gathered to the given
	 * output stream.
	 */
	virtual void dump(std::ostream& os) const = 0;

};
//end class IHandlerRegistry

}//end namespace xolotlPerf

#endif
