#ifndef STDHANDLERREGISTRY_H
#define STDHANDLERREGISTRY_H

#include <map>
#include <memory>
#include "xolotlPerf/IHandlerRegistry.h"


namespace xolotlPerf {

/**
 * Base class for for building performance data collection objects that
 * collect data (as opposed to low-overhead stubs).
 */
class StdHandlerRegistry : public IHandlerRegistry
{
protected:
    /**
     * Collection of the Timers we have created, keyed by name.
     */
    std::map<std::string, std::shared_ptr<ITimer> > allTimers;

    /**
     * Collection of the EventCounters we created, keyed by name.
     */
    std::map<std::string, std::shared_ptr<IEventCounter> > allEventCounters;

    /**
     * Collection of the HWCounterSets we created, keyed by name.
     */
    std::map<std::string, std::shared_ptr<IHardwareCounter> > allHWCounterSets;

    
public:

    /**
     * Construct a StdHandlerRegistry.
     */
    StdHandlerRegistry(void);


    /**
     * Destroy a StdHandlerRegistry.
     */
    virtual ~StdHandlerRegistry(void);


	/**
     * Look up and return a named counter in the current scope.
     * Create the counter if it does not already exist.
     *
     * @param name The object's name.
     * @return The object with the given name.
	 */
	virtual std::shared_ptr<IEventCounter> getEventCounter( std::string name);
};

} // namespace xolotlPerf

#endif // STDHANDLERREGISTRY_H
