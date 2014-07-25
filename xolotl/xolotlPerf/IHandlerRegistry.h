#ifndef IHANDLERREGISTRY_H
#define IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "ITimer.h"
#include "IEventCounter.h"
#include "IHardwareCounter.h"


namespace xolotlPerf {

/**
 * Factory for building performance data collection objects, such
 * as timers and counters.
 */
class IHandlerRegistry {

public:

    /// Possible types of performance handler registries.
    enum RegistryType
    {
        dummy,      //< Use stub classes that do not collect any performance data
        std,        //< Use the best available API.
        os,         //< Use operating system/runtime API.
        papi,       //< Use PAPI to collect performance data.
    };

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
                        const IHardwareCounter::SpecType& ctrSpec ) = 0;

    /**
     * Report statistics about any performance data collected 
     * to the given stream.
     */
    virtual void reportStatistics(std::ostream& os) const = 0;


};

} //end namespace xolotlPerf

#endif
