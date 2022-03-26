#pragma once

#include <memory>
#include <string>

#include <xolotl/perf/IEventCounter.h>
#include <xolotl/perf/IHardwareCounter.h>
#include <xolotl/perf/ITimer.h>
#include <xolotl/perf/PerfObjStatistics.h>

namespace xolotl
{
namespace perf
{
/**
 * Interface for performance handlers
 */
class IPerfHandler
{
public:
	/**
	 * The destructor
	 */
	virtual ~IPerfHandler()
	{
	}

	/**
	 * This operation returns the ITimer specified by the parameter.
	 */
	virtual std::shared_ptr<ITimer>
	getTimer(const std::string& name) = 0;

	/**
	 * This operation returns the IEventCounter specified by the parameter.
	 */
	virtual std::shared_ptr<IEventCounter>
	getEventCounter(const std::string& name) = 0;

	/**
	 * This operation returns the specified IHardwareCounter.
	 */
	virtual std::shared_ptr<IHardwareCounter>
	getHardwareCounter(
		const std::string& name, const IHardwareCounter::SpecType& ctrSpec) = 0;

	/**
	 * Collect statistics about any performance data collected by
	 * processes of the program.
	 *
	 * @param timerStats Map of timer statistics, keyed by timer name.
	 * @param counterStats Map of counter statistics, keyed by counter name.
	 * @param hwCounterStats Map of hardware counter statistics, keyed by
	 * IHardwareCounter name + ':' + hardware counter name.
	 */
	virtual void
	collectStatistics(PerfObjStatsMap<ITimer::ValType>& timerStats,
		PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats) = 0;

	/**
	 * Report performance data statistics to the given stream.
	 *
	 * @param os Stream on which to output statistics.
	 * @param timerStats Map of timer statistics, keyed by timer name.
	 * @param counterStats Map of counter statistics, keyed by counter name.
	 * @param hwCounterStats Map of hardware counter statistics, keyed by
	 * IHardwareCounter name + ':' + hardware counter name.
	 */
	virtual void
	reportStatistics(std::ostream& os,
		const PerfObjStatsMap<ITimer::ValType>& timerStats,
		const PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		const PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats)
		const = 0;
};

void
loadPerfHandlers();
} // namespace perf
} // namespace xolotl
