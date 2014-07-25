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
private:
    template<class T>
    struct PerfObjStatistics
    {
        std::string name;
        unsigned int processCount;
        T min;
        T max;
        double average;
        double stdev;

        PerfObjStatistics( std::string _name )
          : name( _name ),
            processCount( 0 ),
            min( 0 ),
            max( 0 ),
            average( 0 ),
            stdev( 0 )
        { }

        PerfObjStatistics( const PerfObjStatistics& obj )
          : name( obj.name ),
            processCount( obj.processCount ),
            min( obj.min ),
            max( obj.max ),
            average( obj.average ),
            stdev( obj.stdev )
        { }

        PerfObjStatistics& operator=( const PerfObjStatistics& obj )
        {
            if( &obj != this )
            {
                name = obj.name;
                processCount = obj.processCount;
                min = obj.min;
                max = obj.max;
                average = obj.average;
                stdev = obj.stdev;
            }
            return *this;
        }

        void outputTo(std::ostream& os ) const;
    };

    void AggregateTimerStatistics( int myRank, 
            std::map<std::string, PerfObjStatistics<ITimer::ValType> >& timerStats ) const;

    void CollectTimerNames( int myRank,
            std::map<std::string, PerfObjStatistics<ITimer::ValType> >& timerStats ) const;

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


    /**
     * Report statistics about any collected performance data 
     * to the given stream.
     */
    virtual void reportStatistics(std::ostream& os) const;
};

} // namespace xolotlPerf

#endif // STDHANDLERREGISTRY_H
