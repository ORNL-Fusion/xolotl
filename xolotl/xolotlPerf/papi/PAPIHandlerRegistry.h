#ifndef PAPIHANDLERREGISTRY_H
#define PAPIHANDLERREGISTRY_H

#include "xolotlPerf/standard/StdHandlerRegistry.h"


namespace xolotlPerf {

/**
 * Factory for building performance data collection objects that 
 * use the Performance API (PAPI) library for access to timers and
 * hardware counters.
 */
class PAPIHandlerRegistry : public StdHandlerRegistry
{
public:

    /// Construct a handler registry.
    PAPIHandlerRegistry( void );

    /// Destroy the handler registry.
    virtual ~PAPIHandlerRegistry( void ) { }


    /**
     * Look up and return a named timer.
     * Create the timer if it does not already exist.
     *
     * @param name The object's name.
     * @return The object with the given name.
     */
	virtual std::shared_ptr<ITimer> getTimer(const std::string& name);


    /**
     * Look up and return a hardware counter set in the current scope.
     * Create the event counter set if it does not already exist.
     *
     * @param name The object's name.
     * @return The object with the given name.
     */
	virtual std::shared_ptr<IHardwareCounter> getHardwareCounter( const std::string& name,
                        const IHardwareCounter::SpecType& ctrSpec );
};

} // namespace xolotlPerf

#endif // PAPIHandlerRegistry
