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
    /**
     * Statistics (min,max,avg,stdev) for a performance metric we collected
     * during a program run.
     */
    template<class T>
    struct PerfObjStatistics
    {
        std::string name;           //< name of the metric
        unsigned int processCount;  //< number of procesess that collected data for this metric
        T min;                      //< min value across all processes that collected data for the metric
        T max;                      //< max value across all processes that collected data for the metric
        double average;             //< average value across all processes that collected data for the metric
        double stdev;               //< standard deviation of values across all processes that collected data for the metric


        /**
         * Construct a PerfObjStatistics struct with default values.
         * @param _name The metric name.
         */
        PerfObjStatistics( std::string _name )
          : name( _name ),
            processCount( 0 ),
            min( 0 ),
            max( 0 ),
            average( 0 ),
            stdev( 0 )
        { }

        /**
         * Construct a PerfObjStatistics struct as a copy of another.
         * @param obj The object to be copied.
         */
        PerfObjStatistics( const PerfObjStatistics& obj )
          : name( obj.name ),
            processCount( obj.processCount ),
            min( obj.min ),
            max( obj.max ),
            average( obj.average ),
            stdev( obj.stdev )
        { }

        /**
         * Replace my own values to be a copy of another PerfObjStatistics.
         * @param obj The object to be copied.
         */
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

        /**
         * Output our name and statistics to the given output stream.
         * @param os The output stream on which we will write our statistics.
         */
        void outputTo(std::ostream& os ) const;
    };


    /**
     * Collect performance data from all program processes.
     * In the process with rank 0, compute statistics for each
     * performance metric and populate the stats map with those
     * statistics.
     *
     * @param myRank This process' MPI rank.
     * @param myObjs My map of performance data objects of type T (where T
     * is something like ITimer).
     * @param stats A map of statistics structs to be populated with
     * statistics for all performance data objects of type T across all
     * processes in the program.  Only meaningful in process with MPI rank 0.
     */
    template<typename T, typename V>
    void AggregateStatistics( int myRank,
            const std::map<std::string, std::shared_ptr<T> >& myObjs,
            std::map<std::string, PerfObjStatistics<V> >& stats ) const;

    /**
     * Retrieve the value of the named object from the given collection of
     * objects.
     * Note that for hardware counters, the collection should be a collection
     * of IHardwareCounter objects but the name is of the form
     * objName:ctrName where objName is the name associated with an
     * IHardwareCounter and ctrName is a string provided
     * from IHardwareCounter::getCounterName for some 
     * IHardwareCounter::CounterSpec value.
     *
     * @param myObjs A collection of performance data collection objects.
     * @param objName The name of one of the objects in the collection.
     * @return The value of the named object.
     */
    template<typename T, typename V>
    std::pair<bool, V> GetObjValue( const std::map<std::string, std::shared_ptr<T> >& myObjs,
        std::string objName ) const;


    /**
     * Build a collection of names of the performance objects of type T
     * from within my own process.  Type T is something like ITimer.
     *
     * @param myObjs A map of performance data collection objects.
     * @param objNames A vector of performance metric names defined by 
     * the objects in myObjs.
     */
    template<typename T>
    void CollectMyObjectNames( const std::map<std::string,
                                    std::shared_ptr<T> >& myObjs,
                                    std::vector<std::string>& objNames ) const;


    /**
     * Build a collection of all known performance data metric names across 
     * all processes in the program.
     *
     * @param myRank My MPI rank.
     * @param myObjs A collection of the performance data collection objects from my own process.
     * @param stats A map of partially constructed PerfObjStatistics, keyed 
     * by performance metric name.
     * There is one PerfObjStatistics item in the map for each performance
     * metric known across all processes of the program.
     * This map will only be populated within the process with MPI rank 0.
     */
    template<typename T, typename V>
    void CollectAllObjectNames( int myRank,
                    const std::map<std::string, std::shared_ptr<T> >& myObjs,
                    std::map<std::string, PerfObjStatistics<V> >& stats ) const;

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
