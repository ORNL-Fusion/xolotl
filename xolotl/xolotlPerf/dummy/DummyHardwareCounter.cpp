#include "DummyHardwareCounter.h"

using namespace xolotlPerf;


std::vector<double> DummyHardwareCounter::getValues(void) const {
    return std::vector<double>();
}

std::vector<std::string> DummyHardwareCounter::getHardwareQuantities() const {

	return std::vector<std::string>();
}
