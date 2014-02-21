#include "HardwareCounter.h"
#include "HardwareQuantities.h"

using namespace xolotlPerf;

//HardwareCounter::HardwareCounter() :
//		name(""), quantities(hquantities), values(0) {
//}

//HardwareCounter::HardwareCounter() {
//	name = "";
//	values = std::vector<long long>(quantities.size(), 0);
//}

//HardwareCounter::HardwareCounter(){
//
//}

HardwareCounter::HardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities) :
		name(aname), quantities(hquantities), values(std::vector<long long>(quantities.size(), 0)) {
}

HardwareCounter::HardwareCounter(const HardwareCounter &other) :
		name(other.name), quantities(other.quantities), values(other.values) {
}

HardwareCounter::~HardwareCounter() {

}

//std::shared_ptr<HardwareCounter> HardwareCounter::clone() {
//	std::shared_ptr<HardwareCounter> HardwareCounter(new HardwareCounter(*this));
//	return HardwareCounter;
//}

std::vector<long long> HardwareCounter::getValues() const {
	// By default the values are set to zero.
//	int valuesLength = quantities.size();
//	std::vector<long long> values = std::vector<long long>(valuesLength, 0);

	//GPTLget_value to get PAPI counters

	return values;
}

/**
 * This operation returns the name.
 * @return the name
 */
const std::string HardwareCounter::getName() const {
	return name;
}
