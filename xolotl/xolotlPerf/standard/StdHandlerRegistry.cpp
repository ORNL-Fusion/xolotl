#include <iostream>
#include <cassert>
#include "mpi.h"
#include <unistd.h>
#include <float.h>
#include <math.h>
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


template<class T>
void
StdHandlerRegistry::PerfObjStatistics<T>::outputTo(std::ostream& os) const
{
    os << "  " << name << '\n'
        << "    " << "process_count: " << processCount << '\n'
        << "    " << "min: " << min << '\n'
        << "    " << "max: " << max << '\n'
        << "    " << "average: " << average << '\n'
        << "    " << "stdev: " << stdev << '\n'
        << std::endl;
}


void
StdHandlerRegistry::CollectTimerNames( int myRank,
    std::map<std::string, PerfObjStatistics<ITimer::ValType> >& timerStats ) const
{
    // Determine amount of space required for timer names
    unsigned int nBytes = 0;
    for( auto timerIter = allTimers.begin(); timerIter != allTimers.end(); ++timerIter )
    {
        // Add enough space for the timer name plus a NUL terminating character.
        nBytes += (timerIter->first.length() + 1);
    }
    
    // Let root know how much space it needs to collect all timer names
    unsigned int totalNumBytes = 0;
    MPI_Reduce( &nBytes, &totalNumBytes, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

    // Marshal all our timer names.
    char* myTimerNames = new char[nBytes];
    char* pTimerName = myTimerNames;
    for( auto timerIter = allTimers.begin(); timerIter != allTimers.end(); ++timerIter )
    {
        strcpy( pTimerName, timerIter->first.c_str() );
        pTimerName += (timerIter->first.length() + 1);   // skip the NUL terminator
    }
    assert( pTimerName == (myTimerNames + nBytes) );

    // Provide all timer names to root.
    int cwSize;
    MPI_Comm_size( MPI_COMM_WORLD, &cwSize );
    char* allTimerNames = (myRank == 0) ? new char[totalNumBytes] : NULL;
    int* allTimerNameCounts = (myRank == 0) ? new int[cwSize] : NULL;
    int* allTimerNameDispls = (myRank == 0) ? new int[cwSize] : NULL;

    MPI_Gather( &nBytes, 1, MPI_INT, allTimerNameCounts, 1, MPI_INT, 0, MPI_COMM_WORLD );

    if( myRank == 0 )
    {
        allTimerNameDispls[0] = 0;
        for( unsigned int i = 1; i < cwSize; ++i )
        {
            allTimerNameDispls[i] = allTimerNameDispls[i-1] + allTimerNameCounts[i-1];
        }
    }

    MPI_Gatherv( myTimerNames, 
                    nBytes, 
                    MPI_CHAR,
                    allTimerNames,
                    allTimerNameCounts,
                    allTimerNameDispls,
                    MPI_CHAR,
                    0,
                    MPI_COMM_WORLD );

    if( myRank == 0 )
    {
        // Process the gathered timer names to determine the 
        // set of all known timer names.
        pTimerName = allTimerNames;
        while( pTimerName < (allTimerNames + totalNumBytes) )
        {
            auto iter = timerStats.find(pTimerName);
            if( iter == timerStats.end() )
            {
                // This is a timer name we have not seen before.
                // Add it to the timer statistics map.
                timerStats.insert( std::pair<std::string, PerfObjStatistics<ITimer::ValType> >(pTimerName, PerfObjStatistics<ITimer::ValType>(pTimerName) ) );
            }
            
            // Advance to next timer name
            pTimerName += (strlen(pTimerName) + 1);
        }
        assert( pTimerName == allTimerNames + totalNumBytes );
    }

    // clean up
    delete[] myTimerNames;
    delete[] allTimerNames;
    delete[] allTimerNameCounts;
    delete[] allTimerNameDispls;
}


void
StdHandlerRegistry::AggregateTimerStatistics( int myRank,
    std::map<std::string, PerfObjStatistics<ITimer::ValType> >& timerStats ) const
{
    // Determine the set of timer names known across all processes.
    // Since some processes may define a timer that others don't, we
    // have to form the union across all processes.
    // Unfortunately, because the strings are of different lengths,
    // we have a more difficult marshal/unmarshal problem than we'd like.
    CollectTimerNames( myRank, timerStats );

    // Now collect statistics for each timer the program defined.
    int nTimers;
    if( myRank == 0 )
    {
        nTimers = timerStats.size();
    }
    MPI_Bcast( &nTimers, 1, MPI_INT, 0, MPI_COMM_WORLD );
    assert( nTimers >= 0 );

    auto tsiter = timerStats.begin();
    for( unsigned int timerIdx = 0; timerIdx < nTimers; ++timerIdx )
    {
        // broadcast the current timer's name
        int nameLen = (myRank == 0) ? tsiter->second.name.length() : -1;
        MPI_Bcast( &nameLen, 1, MPI_INT, 0, MPI_COMM_WORLD );
        // we can safely cast away const on the tsiter data string because
        // the only process that accesses that string is rank 0,
        // and it only reads the data.
        char* timerName = (myRank == 0) ? const_cast<char*>(tsiter->second.name.c_str()) : new char[nameLen+1];
        MPI_Bcast( timerName, nameLen+1, MPI_CHAR, 0, MPI_COMM_WORLD );

        // do we know about the current timer?
        auto currTimerIter = allTimers.find(timerName);
        int knowTimer = (currTimerIter != allTimers.end()) ? 1 : 0;

        // collect count of processes knowing the current timer
        unsigned int* pcount = (myRank == 0) ? &(tsiter->second.processCount) : NULL;
        MPI_Reduce(&knowTimer, pcount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        // collect min value of current timer        
        ITimer::ValType* pMinVal = (myRank == 0) ? &(tsiter->second.min) : NULL;
        ITimer::ValType myVal = knowTimer ? currTimerIter->second->getValue() : DBL_MAX;
        MPI_Reduce(&myVal, pMinVal, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        // collect max value of current timer
        ITimer::ValType* pMaxVal = (myRank == 0) ? &(tsiter->second.max) : NULL;
        myVal = knowTimer ? currTimerIter->second->getValue() : 0.0;
        MPI_Reduce(&myVal, pMaxVal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // collect sum of current timer's values (for computing avg and stdev)
        double valSum;
        // use the same myVal as for max: actual value if known, 0 otherwise
        MPI_Reduce(&myVal, &valSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if( myRank == 0 )
        {
            tsiter->second.average = valSum / tsiter->second.processCount;
        }

        // collect sum of squares of current timer's values (for stdev)
        double valSquaredSum;
        double myValSquared = myVal*myVal;
        MPI_Reduce(&myValSquared, &valSquaredSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if( myRank == 0 )
        {
            tsiter->second.stdev = sqrt( (valSquaredSum / tsiter->second.processCount) - (tsiter->second.average * tsiter->second.average) );
        }

        // clean up
        if( myRank != 0 )
        {
            delete[] timerName;
        }

        // advance to next timer
        if( myRank == 0 )
        {
            ++tsiter;
        }
    }
}




void
StdHandlerRegistry::reportStatistics(std::ostream& os) const
{
    int myRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

    // Compute statistics about performance data collected by all processes.
    std::map<std::string, PerfObjStatistics<ITimer::ValType> > timerStats;
    std::map<std::string, PerfObjStatistics<IEventCounter::ValType> > counterStats;
    std::map<std::string, PerfObjStatistics<IHardwareCounter::CounterType> > hwCounterStats;
    AggregateTimerStatistics( myRank, timerStats );


    // If I am the rank 0 process, output statistics on the given stream.
    if( myRank == 0 )
    {
        os << "\nTimers:\n";
        for( auto iter = timerStats.begin(); iter != timerStats.end(); ++iter )
        {
            iter->second.outputTo(os);
        }

        os << "\nCounters:\n";
        for( auto iter = counterStats.begin(); iter != counterStats.end(); ++iter )
        {
            iter->second.outputTo(os);
        }

        os << "\nHardwareCounters:\n";
        for( auto iter = hwCounterStats.begin(); iter != hwCounterStats.end(); ++iter )
        {
            iter->second.outputTo(os);
        }
    }
}


} // namespace xolotlPerf

