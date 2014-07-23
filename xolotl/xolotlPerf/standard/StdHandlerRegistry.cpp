#include <iostream>
#include "xolotlPerf/standard/StdHandlerRegistry.h"
#include "xolotlPerf/standard/EventCounter.h"


namespace xolotlPerf
{

StdHandlerRegistry::StdHandlerRegistry(void)
{
    // nothing else to do
}


StdHandlerRegistry::~StdHandlerRegistry(void)
{
    std::cerr << "Destroying the StdHandlerRegistry" << std::endl;

    // Release the objects we have been tracking.
    // Because we use shared_ptrs for these objects,
    // we do not need to explicitly delete the objects themselves.
    allTimers.clear();
    allEventCounters.clear();
    allHWCounterSets.clear();
}


// We can create the EventCounters, since they don't depend on
// more specialized functionality from any of our subclasses.
std::shared_ptr<IEventCounter>
StdHandlerRegistry::getEventCounter(std::string name)
{
    // TODO - associate the object we create with the current region
    std::shared_ptr<IEventCounter> ret;

    // Check if we have already created an event counter with this name.
    auto iter = allEventCounters.find(name);
    if( iter != allEventCounters.end() )
    {
        // We have already created an event counter with this name.
        // Return name.
        ret = iter->second;
    }
    else
    {
        // We have not yet created an event counter with this name.
        // Build one and keep track of it.
        ret = std::make_shared<EventCounter>(name);
        allEventCounters[name] = ret;
    }
    return ret;
}


} // namespace xolotlPerf

