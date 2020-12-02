#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/perf/dummy/DummyEventCounter.h>
#include <xolotl/perf/dummy/DummyHandler.h>
#include <xolotl/perf/dummy/DummyHardwareCounter.h>
#include <xolotl/perf/dummy/DummyTimer.h>

namespace xolotl
{
namespace perf
{
namespace dummy
{
namespace detail
{
auto dummyHandlerRegistrations =
	::xolotl::factory::perf::PerfHandlerFactory::RegistrationCollection<
		DummyHandler>({"dummy"});
}

DummyHandler::DummyHandler(const options::IOptions&)
{
}

// Obtain a Timer by name.
std::shared_ptr<ITimer>
DummyHandler::getTimer(const std::string& name)
{
	return std::make_shared<DummyTimer>();
}

// Obtain an EventCounter by name.
std::shared_ptr<IEventCounter>
DummyHandler::getEventCounter(const std::string& name)
{
	return std::make_shared<DummyEventCounter>();
}

// Obtain a HardwareCounter object by name and by the
// counter data it collects.
std::shared_ptr<IHardwareCounter>
DummyHandler::getHardwareCounter(
	const std::string& name, const IHardwareCounter::SpecType& ctrSpec)
{
	return std::make_shared<DummyHardwareCounter>(ctrSpec);
}

void
DummyHandler::collectStatistics(PerfObjStatsMap<ITimer::ValType>&,
	PerfObjStatsMap<IEventCounter::ValType>&,
	PerfObjStatsMap<IHardwareCounter::CounterType>&)
{
	// do nothing
	return;
}

void
DummyHandler::reportStatistics(std::ostream&,
	const PerfObjStatsMap<ITimer::ValType>&,
	const PerfObjStatsMap<IEventCounter::ValType>&,
	const PerfObjStatsMap<IHardwareCounter::CounterType>&) const
{
	// do nothing
	return;
}
} // namespace dummy
} // namespace perf
} // namespace xolotl
