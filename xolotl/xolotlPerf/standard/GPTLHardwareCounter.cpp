#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include "/home/cxj/Libraries/GPTL/GPTL/gptl-5.0/include/gptl.h"
#include "GPTLHardwareCounter.h"

using namespace xolotlPerf;

GPTLHardwareCounter::~GPTLHardwareCounter() 
{
    //TODO Auto-generated method stub
}

//long long GPTLHardwareCounter::getValue() const {
//
//	double papival = 0.0;
//
//	// The following documentation was taken directly from gptl.c
//	/*
//	** GPTLget_eventvalue: return PAPI-based event value for a timer. All values will be
//	** returned as doubles, even if the event is not derived.
//	**
//	** Input args:
//	** const char *timername: timer name
//	** const char *eventname: event name (must be currently enabled)
//	** int t: thread number (if < 0, the request is for the current thread)
//	**
//	** Output args:
//	** double *value: current value of the event for this timer
//	*/
////	int gret = GPTLget_eventvalue( name.c_str(), -1, &papival );
//
//	return (long long)papival;
//}
