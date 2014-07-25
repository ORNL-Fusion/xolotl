#ifndef IHARDWARECOUNTER_H
#define IHARDWARECOUNTER_H

#include "mpi.h"
#include <vector>
#include "../xolotlCore/IIdentifiable.h"

using namespace std;

namespace xolotlPerf{

/**
 * A counter that can sample hardware counter values.
 */
class IHardwareCounter : public virtual xolotlCore::IIdentifiable {

public:

    /// The types of hardware counters we can collect
    enum CounterSpec
    {
        Instructions,
        Cycles,
        FPOps,
        FPInstructions,
        L1CacheMisses,
        L2CacheMisses,
        L3CacheMisses,
        BranchMispredictions,
    };

    /// Type of an individual IHardwareCounter value.
    typedef long long CounterType;

    /// Type of the IHardwareCounter's collection of values.
    typedef std::vector<CounterType> ValType;

    /// Type of the IHardwareCounter's configuration
    typedef std::vector<CounterSpec> SpecType;


	/**
	 * Destroy th ecounter set.
	 */
	virtual ~IHardwareCounter() { }


    /**
     * Start collecting performance counter data.
     * Throws a xolotlPerf::runtime_error if already collecting.
     */
    virtual void start(void) = 0;


    /**
     * Stop collecting performance counter data.
     * Throws a xolotlPerf::runtime_error if not collecting.
     */
    virtual void stop(void) = 0;


    /**
     * Retrieve the values of the hardware counters that have been collected.
     * The values are only valid if the counter set is not currently counting.
     * @return The current values for our configured hardware counters.
     */
    virtual const ValType& getValues(void) const = 0;


    /**
     * Retrieve the hardware counters we monitor.
     * @return The collection of hardware counters we are configured to monitor.
     */
    virtual const SpecType& getSpecification(void) const = 0;


    /**
     * Retrieve the name of the given hardware counter.
     * @return The name of the given hardware counter.
     */
    virtual std::string getCounterName(IHardwareCounter::CounterSpec cs) const = 0;
};

}  //end namespace xolotlPerf

#endif
