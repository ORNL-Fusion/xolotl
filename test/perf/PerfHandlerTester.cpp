#define BOOST_TEST_MODULE Regression

#include <math.h>
#include <unistd.h>

#include <iostream>
#include <sstream>
#include <string>

#include <boost/test/included/unit_test.hpp>

#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/perf/PerfHandlerRegistry.h>
#include <xolotl/perf/PerfObjStatistics.h>
#include <xolotl/test/MPIFixture.h>

using namespace xolotl;

/**
 * Test suite for HandlerRegistry classes (mainly StdHandler).
 */
BOOST_AUTO_TEST_SUITE(StdHandler_testSuite)

#if BOOST_VERSION >= 105900
// In Boost 1.59, the semicolon at the end of the definition of
// BOOST_GLOBAL_FIXTURE is removed
BOOST_GLOBAL_FIXTURE(MPIFixture);
#else
// With earlier Boost versions, naively adding a semicolon to our code will
// generate compiler warnings about redundant semicolons
BOOST_GLOBAL_FIXTURE(MPIFixture)
#endif

BOOST_AUTO_TEST_CASE(createDummyHandlerReg)
{
	unsigned int nGoodInits = 0;

	try {
		perf::PerfHandlerRegistry::set(
			factory::perf::PerfHandlerFactory::get().generate("dummy"));
		nGoodInits++;

		auto reg = perf::PerfHandlerRegistry::get();
		if (reg) {
			nGoodInits++;
		}

		BOOST_TEST_MESSAGE("Dummy handler registry created successfully.");
	}
	catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
			"DummyHandlerRegistry creation failed: " << e.what());
	}

	BOOST_REQUIRE_EQUAL(nGoodInits, 2U);
}

BOOST_AUTO_TEST_CASE(createOSHandlerReg)
{
	unsigned int nGoodInits = 0;

	try {
		perf::PerfHandlerRegistry::set(
			factory::perf::PerfHandlerFactory::get().generate("os"));
		nGoodInits++;

		auto reg = perf::PerfHandlerRegistry::get();
		if (reg) {
			nGoodInits++;
		}

		BOOST_TEST_MESSAGE("OS handler registry created successfully.");
	}
	catch (const std::exception& e) {
		BOOST_TEST_MESSAGE("OSHandlerRegistry creation failed: " << e.what());
	}

	BOOST_REQUIRE_EQUAL(nGoodInits, 2U);
}

BOOST_AUTO_TEST_CASE(aggregateStats)
{
	try {
		perf::PerfHandlerRegistry::set(
			factory::perf::PerfHandlerFactory::get().generate("os"));
		auto reg = perf::PerfHandlerRegistry::get();

		std::shared_ptr<perf::IEventCounter> ctr =
			reg->getEventCounter("testCounter");
		if (!ctr) {
			throw std::runtime_error("Failed to create EventCounter");
		}

		std::shared_ptr<perf::ITimer> timer = reg->getTimer("testTimer");
		if (!timer) {
			throw std::runtime_error("Failed to create Timer");
		}

		int cwSize = -1;
		int cwRank = -1;
		MPI_Comm_size(MPI_COMM_WORLD, &cwSize);
		MPI_Comm_rank(MPI_COMM_WORLD, &cwRank);

		// simulate the timing of some event.
		BOOST_TEST_MESSAGE("Simulating timed event...");
		const unsigned int nTimedSeconds = 5;
		timer->start();
		sleep(5);
		timer->stop();
		BOOST_TEST_MESSAGE("done.");

		// simulate events that occurred
		BOOST_TEST_MESSAGE("Simulating counting events...");
		const unsigned int nEvents = cwRank;
		for (unsigned int i = 0; i < nEvents; ++i) {
			ctr->increment();
		}
		BOOST_TEST_MESSAGE("done.");

		// compute statistics about the program's event counts
		perf::PerfObjStatsMap<perf::ITimer::ValType> timerStats;
		perf::PerfObjStatsMap<perf::IEventCounter::ValType> ctrStats;
		perf::PerfObjStatsMap<perf::IHardwareCounter::CounterType> hwCtrStats;
		reg->collectStatistics(timerStats, ctrStats, hwCtrStats);

		// Verify the statistics collected.
		// Only rank 0 does the verification.
		if (cwRank == 0) {
			// First check times.  Should be very close to the nTimedSeconds
			// with little spread.
			BOOST_REQUIRE_EQUAL(timerStats.size(), 1U);
			perf::PerfObjStatistics<perf::ITimer::ValType>& timerStatsObj =
				timerStats.begin()->second;

			BOOST_TEST_MESSAGE("timer name: " << timerStatsObj.name);
			BOOST_TEST_MESSAGE(
				"timer process count: " << timerStatsObj.processCount);
			BOOST_TEST_MESSAGE("timer average: " << timerStatsObj.average);
			BOOST_TEST_MESSAGE("timer min: " << timerStatsObj.min);
			BOOST_TEST_MESSAGE("timer max: " << timerStatsObj.max);
			BOOST_TEST_MESSAGE("timer stdev: " << timerStatsObj.stdev);

			BOOST_REQUIRE_EQUAL(timerStatsObj.name, "testTimer");
			BOOST_REQUIRE_EQUAL(
				timerStatsObj.processCount, (unsigned int)cwSize);
			BOOST_REQUIRE_CLOSE(timerStatsObj.average, nTimedSeconds, 0.03);
			BOOST_REQUIRE_CLOSE(timerStatsObj.min, nTimedSeconds, 0.03);
			BOOST_REQUIRE_CLOSE(timerStatsObj.max, nTimedSeconds, 0.03);
			BOOST_REQUIRE_SMALL(timerStatsObj.stdev, 0.001);

			// Next check counts.  We can compute exactly what the values
			// should be.
			double countSum = 0;
			double squaredCountSum = 0;
			for (int i = 0; i < cwSize; ++i) {
				countSum += i;
				squaredCountSum += (i * i);
			}
			double expAverage = (countSum / cwSize);
			unsigned int expMin = 0;
			unsigned int expMax = cwSize - 1;
			double expStdev =
				sqrt((squaredCountSum / cwSize) - expAverage * expAverage);

			BOOST_REQUIRE_EQUAL(ctrStats.size(), 1U);
			perf::PerfObjStatistics<perf::IEventCounter::ValType>& ctrStatsObj =
				ctrStats.begin()->second;

			BOOST_TEST_MESSAGE("ctr name: " << ctrStatsObj.name);
			BOOST_TEST_MESSAGE(
				"ctr process count: " << ctrStatsObj.processCount);
			BOOST_TEST_MESSAGE("ctr average: " << ctrStatsObj.average);
			BOOST_TEST_MESSAGE("ctr min: " << ctrStatsObj.min);
			BOOST_TEST_MESSAGE("ctr max: " << ctrStatsObj.max);
			BOOST_TEST_MESSAGE("ctr stdev: " << ctrStatsObj.stdev);

			BOOST_REQUIRE_EQUAL(ctrStatsObj.name, "testCounter");
			BOOST_REQUIRE_EQUAL(ctrStatsObj.processCount, (unsigned int)cwSize);
			BOOST_REQUIRE_EQUAL(ctrStatsObj.min, expMin);
			BOOST_REQUIRE_EQUAL(ctrStatsObj.max, expMax);
			BOOST_REQUIRE_CLOSE(ctrStatsObj.average, expAverage, 0.01);
			BOOST_REQUIRE_CLOSE(ctrStatsObj.stdev, expStdev, 0.01);
		}
	}
	catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
			"Test of aggregating counter stats failed: " << e.what());
	}
}

BOOST_AUTO_TEST_SUITE_END()
