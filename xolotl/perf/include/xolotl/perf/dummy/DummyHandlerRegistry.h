#ifndef DUMMYHANDLERREGISTRY_H
#define DUMMYHANDLERREGISTRY_H

#include <xolotl/perf/IHandlerRegistry.h>
#include <xolotl/perf/dummy/DummyEventCounter.h> //Dependency Generated Source:DummyHandlerRegistry Target:DummyEventCounter
#include <xolotl/perf/dummy/DummyHardwareCounter.h> //Dependency Generated Source:DummyHandlerRegistry Target:DummyHardwareCounter
#include <xolotl/perf/dummy/DummyTimer.h> //Dependency Generated Source:DummyHandlerRegistry Target:DummyTimer

namespace xolotl
{
namespace perf
{
namespace dummy
{
// Factory for creating timers, event counters, and hardware counter objects
// that are dummies, i.e., they provide the right interface but don't do
// anything.  They are stubs.  This is so that the client code can be
// written to use the performance data collection infrastructure without regard
// to whether performance data collection is active or disabled.
//
class DummyHandlerRegistry : public IHandlerRegistry
{
public:
	DummyHandlerRegistry(void)
	{
	}

	virtual ~DummyHandlerRegistry(void)
	{
	}

	// Obtain a Timer by name.
	virtual std::shared_ptr<ITimer>
	getTimer(const std::string& name);

	// Obtain an EventCounter by name.
	virtual std::shared_ptr<IEventCounter>
	getEventCounter(const std::string& name);

	// Obtain a HardwareCounter object by name and by the
	// counter data it collects.
	virtual std::shared_ptr<IHardwareCounter>
	getHardwareCounter(
		const std::string& name, const IHardwareCounter::SpecType& ctrSpec);

	/**
	 * Collect statistics about any performance data collected by
	 * processes of the program.
	 * This method is a stub.
	 *
	 * @param timerStats Map of timer statistics, keyed by timer name.
	 * @param counterStats Map of counter statistics, keyed by counter name.
	 * @param hwCounterStats Map of hardware counter statistics, keyed by
	 * IHardwareCounter name + ':' + hardware counter name.
	 *
	 */
	virtual void
	collectStatistics(PerfObjStatsMap<ITimer::ValType>& timerStats,
		PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats);

	/**
	 * Report performance data statistics to the given stream.
	 * This method is a stub in this class.
	 *
	 * @param os Stream on which to output statistics.
	 * @param timerStats Map of timer statistics, keyed by timer name.
	 * @param counterStats Map of counter statistics, keyed by counter name.
	 * @param hwCounterStats Map of hardware counter statistics, keyed by
	 * IHardwareCounter name + ':' + hardware counter name.
	 *
	 */
	virtual void
	reportStatistics(std::ostream& os,
		const PerfObjStatsMap<ITimer::ValType>& timerStats,
		const PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		const PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats)
		const;
};

} // end namespace dummy
} // end namespace perf
} // end namespace xolotl

#endif
