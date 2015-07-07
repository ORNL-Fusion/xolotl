#define BOOST_TEST_MODULE Regression

#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <unistd.h>
#include <boost/test/included/unit_test.hpp>
#include "xolotlPerf/xolotlPerf.h"
#include "xolotlPerf/PerfObjStatistics.h"


namespace xperf = xolotlPerf;

// our coordinates in the MPI world
int cwRank = -1;
int cwSize = -1;


/**
 * Test suite for HandlerRegistry classes (mainly StdHandlerRegistry).
 */
BOOST_AUTO_TEST_SUITE (StdHandlerRegistry_testSuite)

struct MPIFixture
{
    MPIFixture( void )
    {
        MPI_Init( &boost::unit_test::framework::master_test_suite().argc, 
                    &boost::unit_test::framework::master_test_suite().argv );

        MPI_Comm_rank( MPI_COMM_WORLD, &cwRank );
        MPI_Comm_size( MPI_COMM_WORLD, &cwSize );
    }

    ~MPIFixture( void )
    {
        MPI_Finalize();
    }
};

BOOST_GLOBAL_FIXTURE( MPIFixture );

BOOST_AUTO_TEST_CASE(createDummyHandlerReg)
{
    unsigned int nGoodInits = 0;

    try
    {
        xperf::initialize( xperf::IHandlerRegistry::dummy );
        nGoodInits++;

        std::shared_ptr<xperf::IHandlerRegistry> reg = xperf::getHandlerRegistry();
        if( reg )
        {
            nGoodInits++;
        }

        BOOST_TEST_MESSAGE("Dummy handler registry created successfully.");
    }
    catch (const std::exception& e)
    {
        BOOST_TEST_MESSAGE("DummyHandlerRegistry creation failed: " << e.what());
    }

    BOOST_REQUIRE_EQUAL( nGoodInits, 2 );
}

BOOST_AUTO_TEST_CASE(createStdHandlerReg)
{
    unsigned int nGoodInits = 0;

    try
    {
        xperf::initialize( xperf::IHandlerRegistry::std );
        nGoodInits++;

        std::shared_ptr<xperf::IHandlerRegistry> reg = xperf::getHandlerRegistry();
        if( reg )
        {
            nGoodInits++;
        }

        BOOST_TEST_MESSAGE("Standard handler registry created successfully.");
    }
    catch (const std::exception& e)
    {
        BOOST_TEST_MESSAGE("StdHandlerRegistry creation failed: " << e.what());
    }

    BOOST_REQUIRE_EQUAL( nGoodInits, 2 );
}

BOOST_AUTO_TEST_CASE(createOSHandlerReg)
{
    unsigned int nGoodInits = 0;

    try
    {
        xperf::initialize( xperf::IHandlerRegistry::os );
        nGoodInits++;

        std::shared_ptr<xperf::IHandlerRegistry> reg = xperf::getHandlerRegistry();
        if( reg )
        {
            nGoodInits++;
        }

        BOOST_TEST_MESSAGE("OS handler registry created successfully.");
    }
    catch (const std::exception& e)
    {
        BOOST_TEST_MESSAGE("OSHandlerRegistry creation failed: " << e.what());
    }

    BOOST_REQUIRE_EQUAL( nGoodInits, 2 );
}


BOOST_AUTO_TEST_CASE(aggregateStats)
{
    try
    {
        xperf::initialize( xperf::IHandlerRegistry::std );
        std::shared_ptr<xperf::IHandlerRegistry> reg = xperf::getHandlerRegistry();

        std::shared_ptr<xperf::IEventCounter> ctr = reg->getEventCounter("testCounter");
        if( !ctr )
        {
            throw std::runtime_error("Failed to create EventCounter");
        }

        std::shared_ptr<xperf::ITimer> timer = reg->getTimer("testTimer");
        if( !timer )
        {
            throw std::runtime_error("Failed to create Timer");
        }

        // simulate the timing of some event.
        BOOST_TEST_MESSAGE("Simulating timed event...");
        const unsigned int nTimedSeconds = 5;
        timer->start();
        sleep( 5 );
        timer->stop();
        BOOST_TEST_MESSAGE("done.");

        // simulate events that occurred
        BOOST_TEST_MESSAGE("Simulating counting events...");
        const unsigned int nEvents = cwRank;
        for( unsigned int i = 0; i < nEvents; ++i )
        {
            ctr->increment();
        }
        BOOST_TEST_MESSAGE("done.");
        
        // compute statistics about the program's event counts
        xperf::PerfObjStatsMap<xperf::ITimer::ValType> timerStats;
        xperf::PerfObjStatsMap<xperf::IEventCounter::ValType> ctrStats;
        xperf::PerfObjStatsMap<xperf::IHardwareCounter::CounterType> hwCtrStats;
        reg->collectStatistics( timerStats, ctrStats, hwCtrStats );

        // Verify the statistics collected.
        // Only rank 0 does the verification.
        if( cwRank == 0 )
        {
            // First check times.  Should be very close to the nTimedSeconds
            // with little spread.
            BOOST_REQUIRE_EQUAL( timerStats.size(), 1 );
            xperf::PerfObjStatistics<xperf::ITimer::ValType>& timerStatsObj = timerStats.begin()->second;

            BOOST_TEST_MESSAGE( "timer name: " << timerStatsObj.name );
            BOOST_TEST_MESSAGE( "timer process count: " << timerStatsObj.processCount );
            BOOST_TEST_MESSAGE( "timer average: " << timerStatsObj.average );
            BOOST_TEST_MESSAGE( "timer min: " << timerStatsObj.min );
            BOOST_TEST_MESSAGE( "timer max: " << timerStatsObj.max );
            BOOST_TEST_MESSAGE( "timer stdev: " << timerStatsObj.stdev );

            BOOST_REQUIRE_EQUAL( timerStatsObj.name, "testTimer" );
            BOOST_REQUIRE_EQUAL( timerStatsObj.processCount, cwSize );
            BOOST_REQUIRE_CLOSE( timerStatsObj.average, nTimedSeconds, 0.03);
            BOOST_REQUIRE_CLOSE( timerStatsObj.min, nTimedSeconds, 0.03);
            BOOST_REQUIRE_CLOSE( timerStatsObj.max, nTimedSeconds, 0.03);
            BOOST_REQUIRE_SMALL( timerStatsObj.stdev, 0.001 );

            // Next check counts.  We can compute exactly what the values
            // should be.
            double countSum = 0;
            double squaredCountSum = 0;
            for( unsigned int i = 0; i < cwSize; ++i )
            {
                countSum += i;
                squaredCountSum += (i*i);
            }
            double expAverage = (countSum / cwSize);
            unsigned int expMin = 0;
            unsigned int expMax = cwSize-1;
            double expStdev = sqrt( (squaredCountSum / cwSize) - expAverage*expAverage);

            BOOST_REQUIRE_EQUAL( ctrStats.size(), 1 );
            xperf::PerfObjStatistics<xperf::IEventCounter::ValType>& ctrStatsObj = ctrStats.begin()->second;

            BOOST_TEST_MESSAGE( "ctr name: " << ctrStatsObj.name );
            BOOST_TEST_MESSAGE( "ctr process count: " << ctrStatsObj.processCount );
            BOOST_TEST_MESSAGE( "ctr average: " << ctrStatsObj.average );
            BOOST_TEST_MESSAGE( "ctr min: " << ctrStatsObj.min );
            BOOST_TEST_MESSAGE( "ctr max: " << ctrStatsObj.max );
            BOOST_TEST_MESSAGE( "ctr stdev: " << ctrStatsObj.stdev );

            BOOST_REQUIRE_EQUAL( ctrStatsObj.name, "testCounter" );
            BOOST_REQUIRE_EQUAL( ctrStatsObj.processCount, cwSize );
            BOOST_REQUIRE_EQUAL( ctrStatsObj.min, expMin );
            BOOST_REQUIRE_EQUAL( ctrStatsObj.max, expMax );
            BOOST_REQUIRE_CLOSE( ctrStatsObj.average, expAverage, 0.01 );
            BOOST_REQUIRE_CLOSE( ctrStatsObj.stdev, expStdev, 0.01 );
        }
    }
    catch (const std::exception& e)
    {
        BOOST_TEST_MESSAGE("Test of aggregating counter stats failed: " << e.what());
    }
}

BOOST_AUTO_TEST_SUITE_END()
