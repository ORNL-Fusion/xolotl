#include <xolotl/perf/os/OSHandlerRegistry.h>
#include <xolotl/perf/os/OSTimer.h>
#include <xolotl/perf/EventCounter.h>
#include <xolotl/perf/dummy/DummyHardwareCounter.h>

namespace xolotlPerf {

std::shared_ptr<ITimer> OSHandlerRegistry::getTimer(const std::string& name) {
	// TODO - associate the object we create with the current region.
	std::shared_ptr<ITimer> ret;

	// check if we have already created a timer with this name
	auto iter = allTimers.find(name);
	if (iter != allTimers.end()) {
		// We have already created a timer with this name.
		// Return it.
		ret = iter->second;
	} else {
		// We have not yet created a timer with this name.
		// Build one, and keep track of it.
		ret = std::make_shared<OSTimer>(name);
		allTimers[name] = ret;
	}
	return ret;
}

std::shared_ptr<IHardwareCounter> OSHandlerRegistry::getHardwareCounter(
		const std::string& name, const IHardwareCounter::SpecType& ctrSpec) {
	// TODO - associate the object we create with the current region
	std::shared_ptr<IHardwareCounter> ret;

	// Check if we have already created a dummy hardware counter set with this name.
	auto iter = allHWCounterSets.find(name);
	if (iter != allHWCounterSets.end()) {
		// We have already created a hw counter set with this name.
		// Return it.
		ret = iter->second;
	} else {
		// We have not yet created a hw counter set with this name.
		// Build one and keep track of it.
		// Note with the OSHandlerRegistry it is always a dummy.
		ret = std::make_shared<DummyHardwareCounter>(name, ctrSpec);
		allHWCounterSets[name] = ret;
	}
	return ret;
}

} // namespace xolotlPerf

