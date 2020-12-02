#include <papi.h>

#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/perf/EventCounter.h>
#include <xolotl/perf/RuntimeError.h>
#include <xolotl/perf/papi/PAPIHandler.h>
#include <xolotl/perf/papi/PAPIHardwareCounter.h>
#include <xolotl/perf/papi/PAPITimer.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
namespace detail
{
auto papiHandlerRegistrations =
	::xolotl::factory::perf::PerfHandlerFactory::RegistrationCollection<
		PAPIHandler>({"papi"});
}

PAPIHandler::PAPIHandler(const options::IOptions& options) :
	PerfHandler(options)
{
	int ret;
	ret = PAPI_library_init(PAPI_VER_CURRENT);
	if (ret != PAPI_VER_CURRENT) {
		throw xolotl::perf::runtime_error(
			"Unable to initialize PAPI library for performance data collection",
			ret);
	}
}

PAPIHandler::~PAPIHandler()
{
}

std::shared_ptr<ITimer>
PAPIHandler::getTimer(const std::string& name)
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
		ret = std::make_shared<PAPITimer>();
		allTimers[name] = ret;
	}
	return ret;
}

std::shared_ptr<IHardwareCounter>
PAPIHandler::getHardwareCounter(
	const std::string& name, const IHardwareCounter::SpecType& ctrSpec)
{
	// TODO - associate the object we create with the current region
	std::shared_ptr<IHardwareCounter> ret;

	// Check if we have already created a hardware counter set with this name.
	auto iter = allHWCounterSets.find(name);
	if (iter != allHWCounterSets.end()) {
		// We have already created a hw counter set with this name.
		// Return it.
		ret = iter->second;
	}
	else {
		// We have not yet created a hw counter set with this name.
		// Build one and keep track of it.
		ret = std::make_shared<PAPIHardwareCounter>(ctrSpec);
		allHWCounterSets[name] = ret;
	}
	return ret;
}

} // namespace papi
} // namespace perf
} // namespace xolotl
