#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H

#include <iostream>
#include <string>
#include <map>
#include "IHandlerRegistry.h"
#include "GPTLTimer.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLTimer
#include "GPTLHardwareCounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLHardwareCounter
#include "EventCounter.h"
#include "HardwareQuantityInfoMap.h"


namespace xolotlPerf {


// Factory for creating timers and hardware counter objects
// objects that use GPTL for collecting performance data,
// and event counters that interoperate with the timers and hardware counter
// objects.
//
// This is named "Standard" because it is expected to be the standard
// performance data collection mechanism used by XOLOTL.
// TODO Perhaps "Default" would be a better name?
class StandardHandlerRegistry : public IHandlerRegistry
{
private:
    // Map of our hardware quantity values to 
    // human readable names and PAPI counter IDs.
    HardwareQuantityInfoMap hwqInfoMap;


    // Collection of the EventCounters we created.
    std::map<std::string,std::shared_ptr<IEventCounter>> allEventCounters;

public:
    // Construct a StandardHandlerRegistry.
    // @param hwq The hardware quantities we will be collecting during this run.
    StandardHandlerRegistry( std::vector<HardwareQuantities> hwq );

    // Clean up a StandardHandlerRegistry.
    virtual ~StandardHandlerRegistry( void );

    // Obtain a Timer by name.
    // @param name The name of the timer.
    // @return A shared pointer to the newly-created timer.
    virtual std::shared_ptr<ITimer> getTimer(std::string name);

    // Create or access an EventCounter by name.
    // @param name The name of the EventCounter.
    // @return A shared pointer to the newly-created event counter.
    virtual std::shared_ptr<IEventCounter> getEventCounter(std::string name);

    // Obtain a HardwareCounter object by name and by the 
    // counter data it collects.
    //
    // NOTE: due to the way that GPTL works, it is incapable of
    // configuring the processor to collect fine-grained data on
    // hardware counters.  Instead, pass all the event counters you
    // wish to collect during the run to the constructor.
    //
    // @param name The name of the hardware counter.
    // @param quantities The hardware counters we will be collecting.
    // @return A shared pointer to the newly-created event counter.
	virtual std::shared_ptr<IHardwareCounter> getHardwareCounter(
			std::string name,
			std::vector<HardwareQuantities> quantities);

	// This operation outputs the information gathered to the given 
    // output stream.
    // @param os The output stream to which performance data is written.
	virtual void dump(std::ostream& os) const;

};  //end class StandardHandlerRegistry

}//end namespace xolotlPerf

#endif
