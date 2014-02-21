#include "Timer.h"
#include <string>

using namespace xolotlPerf;

Timer::Timer() :
		name("") {
}

Timer::Timer(std::string aname) :
		name(aname) {
}

Timer::Timer(const Timer &other) :
		name(other.name) {
}

Timer::~Timer() {

}

//std::shared_ptr<Timer> Timer::clone() {
//	std::shared_ptr<Timer> timer(new Timer(*this));
//	return timer;
//}


/**
 * This operation returns the name.
 * @return the name
 */
const std::string Timer::getName() const {
	return name;
}



