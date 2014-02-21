#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "Timer.h"
#include "GPTLTimer.h"


using namespace xolotlPerf;

//GPTLTimer::GPTLTimer():Timer("") {
//	GPTLinitialize();
//}

GPTLTimer::GPTLTimer(std::string aname) : Timer(aname) {

	value = 0.0;

//~~~~~~~~~~~~~~~~~MOVE SOME WHERE
//	GPTLsetoption (GPTLcpu, 1);
//	GPTLsetoption (GPTLwall, 1);
//	GPTLinitialize();
//	GPTLinit_handle(name.c_str(), &handle);
}

//GPTLTimer::GPTLTimer(const std::string &name) : Timer(name) {
//	GPTLinitialize();
//	GPTLinit_handle(name.c_str(), &handle);
//}

GPTLTimer::~GPTLTimer() {
}

// This operations starts the Timer.
void GPTLTimer::start(){
//	GPTLstart_handle(aname.c_str(), &handle);
	GPTLstart(name.c_str());
}

// This operation stops the Timer.
void GPTLTimer::stop(){
//	GPTLstop_handle(name.c_str(), &handle);
	GPTLstop(name.c_str());
}

const std::string GPTLTimer::getName() const {
	return name;
}

// This operation returns the value of the Timer.
double GPTLTimer::getValue() {

	/*
	** GPTLget_wallclock: return wallclock accumulation for a timer.
	**
	** Input args:
	** const char *timername: timer name
	** int t: thread number (if < 0, the request is for the current thread)
	**
	** Output args:
	** double *value: current wallclock accumulation for the timer
	*/
    int gret = GPTLget_wallclock( name.c_str(), -1, &value );

    return value;
}

// This operation returns the units of the Timer.
long GPTLTimer::getUnits() const {
//	return this->units;
}

