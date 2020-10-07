#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/perf/EventCounter.h>
#include <xolotl/perf/dummy/DummyHardwareCounter.h>
#include <xolotl/perf/os/OSHandler.h>
#include <xolotl/perf/os/OSTimer.h>

namespace xolotl
{
namespace perf
{
namespace os
{
namespace detail
{
auto osHandlerRegistrations =
	::xolotl::factory::perf::PerfHandlerFactory::RegistrationCollection<
		OSHandler>({"os", "std"});
}

OSHandler::OSHandler(const options::IOptions& options) :
	standard::StdHandler(options)
{
}

std::shared_ptr<ITimer>
OSHandler::getTimer(const std::string& name)
{
	// TODO - associate the object we create with the current region.
	std::shared_ptr<ITimer> ret;

	// check if we have already created a timer with this name
	auto iter = allTimers.find(name);
	if (iter != allTimers.end()) {
		// We have already created a timer with this name.
		// Return it.
		ret = iter->second;
	}
	else {
		// We have not yet created a timer with this name.
		// Build one, and keep track of it.
		ret = std::make_shared<OSTimer>(name);
		allTimers[name] = ret;
	}
	return ret;
}
} // namespace os
} // namespace perf
} // namespace xolotl
