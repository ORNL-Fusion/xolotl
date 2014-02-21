#include "DummyTimer.h"
#include <time.h>

using namespace xolotlPerf;


DummyTimer::DummyTimer(std::string aname) : Timer(aname) {

}

DummyTimer::~DummyTimer() {

}

//std::shared_ptr<Timer> DummyTimer::clone() {
//	std::shared_ptr<Timer> timer(new DummyTimer(*this));
//	return timer;
//}

void DummyTimer::start() {

}

void DummyTimer::stop() {

}

const std::string DummyTimer::getName() const {
	return "";
}

double DummyTimer::getValue() {
	return 0;
}

long DummyTimer::getUnits() const {
	return 0;
}
