#include "HardwareQuantities.h"
#include "DummyHardwareCounter.h"

using namespace xolotlPerf;

DummyHardwareCounter::DummyHardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities) :
		HardwareCounter(aname, hquantities) {
}

DummyHardwareCounter::~DummyHardwareCounter() {

}

//std::shared_ptr<HardwareCounter> DummyHardwareCounter::clone() {
//	std::shared_ptr<HardwareCounter> eventCounter(new DummyHardwareCounter(*this));
//	return eventCounter;
//}

std::vector<long long> DummyHardwareCounter::getValues() const {
	return values;
}

/**
 * This operation returns the name.
 * @return the name
 */
const std::string DummyHardwareCounter::getName() const {
	return "";
}

void DummyHardwareCounter::increment(){

}
