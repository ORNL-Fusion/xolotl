#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "GPTLEventCounter.h"

using namespace xolotlPerf;


GPTLEventCounter::GPTLEventCounter(std::string aname) : EventCounter(aname) {

	value = 0;

}

GPTLEventCounter::~GPTLEventCounter() {

}

//GPTLEventCounter::GPTLEventCounter(const GPTLEventCounter &other) :
//		name(other.name), value(other.value) {
//}


int GPTLEventCounter::getValue() {

//	int val = 0;

	// The following documentation was taken directly from gptl.c
	/*
	** GPTLget_eventvalue: return PAPI-based event value for a timer. All values will be
	** returned as doubles, even if the event is not derived.
	**
	** Input args:
	** const char *timername: timer name
	** const char *eventname: event name (must be currently enabled)
	** int t: thread number (if < 0, the request is for the current thread)
	**
	** Output args:
	** double *value: current value of the event for this timer
	*/
//	int gret = GPTLget_eventvalue( name.c_str(), -1, &papival );


	/*
	** GPTLget_count: return number of start/stop calls for a timer.
	**
	** Input args:
	** const char *timername: timer name
	** int t: thread number (if < 0, the request is for the current thread)
	**
	** Output args:
	** int *count: current number of start/stop calls for the timer
	*/
	int gret = GPTLget_count( name.c_str(), -1, &value );

	return value;
}

const std::string GPTLEventCounter::getName() const {
	return name;
}

void GPTLEventCounter::increment(){

	++value;

}
