#include "DummyEventCounter.h"

using namespace xolotlPerf;


DummyEventCounter::DummyEventCounter(std::string aname) : EventCounter(aname) {
}

//DummyEventCounter::DummyEventCounter(const DummyEventCounter &other) :
//		name(other.name), value(other.value) {
//}

DummyEventCounter::~DummyEventCounter() {

}

//std::shared_ptr<EventCounter> DummyEventCounter::clone() {
//	std::shared_ptr<EventCounter> eventCounter(new DummyEventCounter(*this));
//	return eventCounter;
//}

int DummyEventCounter::getValue() {
	return 0;
}

/**
 * This operation returns the name.
 * @return the name
 */
const std::string DummyEventCounter::getName() const {
	return "";
}

void DummyEventCounter::increment(){

}
