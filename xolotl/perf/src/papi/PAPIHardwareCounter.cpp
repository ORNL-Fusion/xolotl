#include <papi.h>

#include <cassert>
#include <iostream>
#include <sstream>

#include <xolotl/perf/RuntimeError.h>
#include <xolotl/perf/papi/PAPIHardwareCounter.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
PAPIHardwareCounter::CounterSpecMap PAPIHardwareCounter::csMap;

PAPIHardwareCounter::PAPIHardwareCounter(
	const IHardwareCounter::SpecType& cset) :
	spec(cset),
	eventSet(PAPI_NULL)
{
	assert(PAPI_is_initialized());

	// Ensure our counter spec map has been initialized.
	if (csMap.empty()) {
		initCounterSpecMap();
	}

	// Build the PAPI event set for our events.
	int err = PAPI_create_eventset(&eventSet);
	if (err != PAPI_OK) {
		throw xolotl::perf::runtime_error(
			"Failed to create PAPI eventset", err);
	}
	for (const auto& cspec : spec) {
		auto miter = csMap.find(cspec);

		// we had better know about the counter spec
		assert(miter != csMap.end());

		const CounterSpecInfo& currCounterSpecInfo = miter->second;
		err = PAPI_add_event(eventSet, currCounterSpecInfo.papiEventID);
		if (err != PAPI_OK) {
			std::string errorMsg = "Failed to add event to PAPI eventset";
			if (PAPI_query_event(currCounterSpecInfo.papiEventID) != PAPI_OK) {
				errorMsg += "\n" + currCounterSpecInfo.name + " (" +
					currCounterSpecInfo.papiName + ") counter unavailable";
			}
			throw xolotl::perf::runtime_error(errorMsg, err);
		}
	}

	// Ensure our value vector is big enough to collect the events we
	// are supposed to be configured for.
	vals.resize(spec.size(), 0);
}

PAPIHardwareCounter::~PAPIHardwareCounter(void)
{
	if (eventSet != PAPI_NULL) {
		PAPI_cleanup_eventset(eventSet);
		PAPI_destroy_eventset(&eventSet);
		eventSet = PAPI_NULL;
	}
}

void
PAPIHardwareCounter::start(void)
{
	int err = PAPI_start(eventSet);
	if (err != PAPI_OK) {
		throw xolotl::perf::runtime_error(
			"Failed to start collecting configured PAPI events", err);
	}
}

void
PAPIHardwareCounter::stop(void)
{
	assert(vals.size() > 0);
	int err = PAPI_stop(eventSet, &(vals.front()));
	if (err != PAPI_OK) {
		throw xolotl::perf::runtime_error(
			"Failed to stop collecting configured PAPI events", err);
	}
}

std::string
PAPIHardwareCounter::getCounterName(IHardwareCounter::CounterSpec cs) const
{
	std::string ret;

	CounterSpecMap::const_iterator iter = csMap.find(cs);
	if (iter != csMap.end()) {
		std::ostringstream ostr;
		ostr << iter->second.name << '(' << iter->second.papiName << ')';
		ret = ostr.str();
	}
	return ret;
}

void
PAPIHardwareCounter::initCounterSpecMap()
{
	csMap.try_emplace(
		Instructions, "Instructions", "PAPI_TOT_INS", PAPI_TOT_INS);

	csMap.try_emplace(Cycles, "Total cycles", "PAPI_TOT_CYC", PAPI_TOT_CYC);

	csMap.try_emplace(
		FPOps, "Floating point operations", "PAPI_FP_OPS", PAPI_FP_OPS);

	csMap.try_emplace(FPInstructions, "Floating point instructions",
		"PAPI_FP_INS", PAPI_FP_INS);

	csMap.try_emplace(
		L1CacheMisses, "L1 cache misses", "PAPI_L1_TCM", PAPI_L1_TCM);

	csMap.try_emplace(
		L2CacheMisses, "L2 cache misses", "PAPI_L2_TCM", PAPI_L2_TCM);

	csMap.try_emplace(
		L3CacheMisses, "L3 cache misses", "PAPI_L3_TCM", PAPI_L3_TCM);

	csMap.try_emplace(BranchMispredictions, "Branch mispredictions",
		"PAPI_BR_MSP", PAPI_BR_MSP);
}

IHardwareCounter&
PAPIHardwareCounter::operator+=(const IHardwareCounter& c)
{
#if READY
#else
	std::cerr << "PAPIHardwareCounter::NIY" << std::endl;
#endif // READY
	return *this;
}

} // namespace papi
} // namespace perf
} // namespace xolotl
