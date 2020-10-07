#pragma once

#include <map>
#include <memory>

#include <xolotl/perf/IPerfHandler.h>
#include <xolotl/perf/PerfObjStatistics.h>

namespace xolotl
{
namespace perf
{
namespace standard
{
/**
 * Base class for for building performance data collection objects that
 * collect data (as opposed to low-overhead stubs).
 */
class StdHandler : public IPerfHandler
{
private:
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
	template <typename T, typename V>
	void
	aggregateStatistics(int myRank,
		const std::map<std::string, std::shared_ptr<T>>& myObjs,
		std::map<std::string, PerfObjStatistics<V>>& stats) const;

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
	template <typename T, typename V>
	std::pair<bool, V>
	getObjValue(const std::map<std::string, std::shared_ptr<T>>& myObjs,
		const std::string& objName) const;

	/**
	 * Build a collection of names of the performance objects of type T
	 * from within my own process.  Type T is something like ITimer.
	 *
	 * @param myObjs A map of performance data collection objects.
	 * @param objNames A vector of performance metric names defined by
	 * the objects in myObjs.
	 */
	template <typename T>
	void
	collectMyObjectNames(
		const std::map<std::string, std::shared_ptr<T>>& myObjs,
		std::vector<std::string>& objNames) const;

	/**
	 * Build a collection of all known performance data metric names across
	 * all processes in the program.
	 *
	 * @param myRank My MPI rank.
	 * @param myObjs A collection of the performance data collection objects
	 * from my own process.
	 * @param stats A map of partially constructed PerfObjStatistics, keyed
	 * by performance metric name.
	 * There is one PerfObjStatistics item in the map for each performance
	 * metric known across all processes of the program.
	 * This map will only be populated within the process with MPI rank 0.
	 */
	template <typename T, typename V>
	void
	collectAllObjectNames(int myRank,
		const std::map<std::string, std::shared_ptr<T>>& myObjs,
		std::map<std::string, PerfObjStatistics<V>>& stats) const;

protected:
	/**
	 * Collection of the Timers we have created, keyed by name.
	 */
	std::map<std::string, std::shared_ptr<ITimer>> allTimers;

	/**
	 * Collection of the EventCounters we created, keyed by name.
	 */
	std::map<std::string, std::shared_ptr<IEventCounter>> allEventCounters;

	/**
	 * Collection of the HWCounterSets we created, keyed by name.
	 */
	std::map<std::string, std::shared_ptr<IHardwareCounter>> allHWCounterSets;

public:
    StdHandler(const options::IOptions&);

	/**
	 * Destroy a StdHandler.
	 */
	virtual ~StdHandler();

	/**
	 * Look up and return a named counter in the current scope.
	 * Create the counter if it does not already exist.
	 *
	 * @param name The object's name.
	 * @return The object with the given name.
	 */
	std::shared_ptr<IEventCounter>
	getEventCounter(const std::string& name) override;

	/**
	 * Look up and return a named counter in the current scope.
	 * Create the counter if it does not already exist.
	 *
	 * @param name The object's name.
	 * @return The object with the given name.
	 */
	std::shared_ptr<IHardwareCounter>
	getHardwareCounter(const std::string& name,
		const IHardwareCounter::SpecType& ctrSpec) override;

	/**
	 * Collect statistics about any performance data collected by
	 * processes of the program.
	 *
	 * @param timerStats Map of timer statistics, keyed by timer name.
	 * @param counterStats Map of counter statistics, keyed by counter name.
	 * @param hwCounterStats Map of hardware counter statistics, keyed by
	 * IHardwareCounter name + ':' + hardware counter name.
	 *
	 */
	void
	collectStatistics(PerfObjStatsMap<ITimer::ValType>& timerStats,
		PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		PerfObjStatsMap<IHardwareCounter::CounterType>& hwStats) override;

	/**
	 * Report performance data statistics to the given stream.
	 *
	 * @param os Stream on which to output statistics.
	 * @param timerStats Map of timer statistics, keyed by timer name.
	 * @param counterStats Map of counter statistics, keyed by counter name.
	 * @param hwCounterStats Map of hardware counter statistics, keyed by
	 * IHardwareCounter name + ':' + hardware counter name.
	 */
	void
	reportStatistics(std::ostream& os,
		const PerfObjStatsMap<ITimer::ValType>& timerStats,
		const PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		const PerfObjStatsMap<IHardwareCounter::CounterType>& hwStats)
		const override;
};
} // namespace standard
} // namespace perf
} // namespace xolotl
