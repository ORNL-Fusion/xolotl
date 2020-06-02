#include <iostream>
#include <sstream>
#include <cassert>
#include <papi.h>
#include <xolotl/perf/papi/PAPIHardwareCounter.h>
#include <xolotl/perf/RuntimeError.h>

namespace xolotl {
namespace perf {
namespace papi {

PAPIHardwareCounter::CounterSpecMap PAPIHardwareCounter::csMap;

PAPIHardwareCounter::PAPIHardwareCounter(const std::string& name,
		const IHardwareCounter::SpecType& cset) :
		core::Identifiable(name), spec(cset), eventSet(PAPI_NULL) {
	assert(PAPI_is_initialized());

	// Ensure our counter spec map has been initialized.
	if (csMap.empty()) {
		InitCounterSpecMap();
	}

	// Build the PAPI event set for our events.
	int err = PAPI_create_eventset(&eventSet);
	if (err != PAPI_OK) {
		throw xolotlPerf::runtime_error("Failed to create PAPI eventset", err);
	}
	for (SpecType::const_iterator iter = spec.begin(); iter != spec.end();
			++iter) {
		CounterSpecMap::const_iterator miter = csMap.find(*iter);

		// we had better know about the counter spec
		assert(miter != csMap.end());

		CounterSpecInfo* currCounterSpecInfo = miter->second;
		err = PAPI_add_event(eventSet, currCounterSpecInfo->papiEventID);
		if (err != PAPI_OK) {
			throw xolotlPerf::runtime_error(
					"Failed to add event to PAPI eventset", err);
		}
	}

	// Ensure our value vector is big enough to collect the events we
	// are supposed to be configured for.
	vals.resize(spec.size(), 0);
}

PAPIHardwareCounter::~PAPIHardwareCounter(void) {
	if (eventSet != PAPI_NULL) {
		PAPI_cleanup_eventset(eventSet);
		PAPI_destroy_eventset(&eventSet);
		eventSet = PAPI_NULL;
	}
}

void PAPIHardwareCounter::start(void) {
	int err = PAPI_start(eventSet);
	if (err != PAPI_OK) {
		throw xolotlPerf::runtime_error(
				"Failed to start collecting configured PAPI events", err);
	}
}

void PAPIHardwareCounter::stop(void) {
	assert(vals.size() > 0);
	int err = PAPI_stop(eventSet, &(vals.front()));
	if (err != PAPI_OK) {
		throw xolotlPerf::runtime_error(
				"Failed to stop collecting configured PAPI events", err);
	}
}

std::string PAPIHardwareCounter::getCounterName(
		IHardwareCounter::CounterSpec cs) const {
	std::string ret;

	CounterSpecMap::const_iterator iter = csMap.find(cs);
	if (iter != csMap.end()) {
		std::ostringstream ostr;
		ostr << iter->second->name << '(' << iter->second->papiName << ')';
		ret = ostr.str();
	}
	return ret;
}

void PAPIHardwareCounter::InitCounterSpecMap(void) {
	csMap[Instructions] = new CounterSpecInfo("Instructions", "PAPI_TOT_INS",
			PAPI_TOT_INS);

	csMap[Cycles] = new CounterSpecInfo("Total cycles", "PAPI_TOT_CYC",
			PAPI_TOT_CYC);

	csMap[FPOps] = new CounterSpecInfo("Floating point operations",
			"PAPI_FP_OPS", PAPI_FP_OPS);

	csMap[FPInstructions] = new CounterSpecInfo("Floating point instructions",
			"PAPI_FP_INS", PAPI_FP_INS);

	csMap[L1CacheMisses] = new CounterSpecInfo("L1 cache misses", "PAPI_L1_TCM",
			PAPI_L1_TCM);

	csMap[L2CacheMisses] = new CounterSpecInfo("L2 cache misses", "PAPI_L2_TCM",
			PAPI_L2_TCM);

	csMap[L3CacheMisses] = new CounterSpecInfo("L3 cache misses", "PAPI_L3_TCM",
			PAPI_L3_TCM);

	csMap[BranchMispredictions] = new CounterSpecInfo("Branch mispredictions",
			"PAPI_BR_MSP", PAPI_BR_MSP);
}

IHardwareCounter&
PAPIHardwareCounter::operator+=(const IHardwareCounter& c) {
#if READY
#else
	std::cerr << "PAPIHardwareCounter::NIY" << std::endl;
#endif // READY
	return *this;
}

} // namespace papi
} // namespace perf
} // namespace xolotl
