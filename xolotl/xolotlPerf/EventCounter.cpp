#include "EventCounter.h"
#include <string>
#include <iostream>

using namespace xolotlPerf;

EventCounter::EventCounter() :
		name("") {
}

EventCounter::EventCounter(std::string aname) :
		name(aname) {
}

EventCounter::EventCounter(const EventCounter &other) :
		name(other.name) {
}

EventCounter::~EventCounter() {

}

//std::shared_ptr<EventCounter> EventCounter::clone() {
//	std::shared_ptr<EventCounter> eventCounter(new EventCounter(*this));
//	return eventCounter;
//}

const std::string EventCounter::getName() const {
	return name;
}
