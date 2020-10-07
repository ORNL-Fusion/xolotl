#pragma once

#include <xolotl/perf/IPerfHandler.h>

namespace xolotl
{
namespace perf
{
namespace dummy
{
class DummyHandler : public IPerfHandler
{
public:
    DummyHandler(const options::IOptions&);

	virtual ~DummyHandler()
	{
	}

	std::shared_ptr<ITimer>
	getTimer(const std::string& name) override;

	std::shared_ptr<IEventCounter>
	getEventCounter(const std::string& name) override;

	std::shared_ptr<IHardwareCounter>
	getHardwareCounter(const std::string& name,
		const IHardwareCounter::SpecType& ctrSpec) override;

	void
	collectStatistics(PerfObjStatsMap<ITimer::ValType>& timerStats,
		PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats)
		override;

	void
	reportStatistics(std::ostream& os,
		const PerfObjStatsMap<ITimer::ValType>& timerStats,
		const PerfObjStatsMap<IEventCounter::ValType>& counterStats,
		const PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats)
		const override;
};
} // namespace dummy
} // namespace perf
} // namespace xolotl
