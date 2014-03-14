#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include "gptl.h"
#include <papi.h>
#include "GPTLHardwareCounter.h"

using namespace xolotlPerf;

// Create the static map of hardware quantities
std::unordered_map<int,hardwareQuantity> GPTLHardwareCounter::AllHardwareQuantitiesMap;

GPTLHardwareCounter::GPTLHardwareCounter(std::string name,
		const std::vector<HardwareQuantities> &counterQuantities)
		: xolotlCore::Identifiable(name), quantities(counterQuantities) {


	static hardwareQuantity l1_cache_miss = {"L1_CACHE_MISS",PAPI_L1_TCM, "PAPI_L1_TCM"};
	static hardwareQuantity l2_cache_miss = {"L2_CACHE_MISS",PAPI_L2_TCM, "PAPI_L2_TCM"};
	static hardwareQuantity l3_cache_miss = {"L3_CACHE_MISS",PAPI_L3_TCM, "PAPI_L3_TCM"};
	static hardwareQuantity branch_mispredictions = {"BRANCH_MISPREDICTIONS",PAPI_BR_MSP, "PAPI_BR_MSP"};
	static hardwareQuantity total_cycles = {"TOTAL_CYCLES",PAPI_TOT_CYC, "PAPI_TOT_CYC"};
	static hardwareQuantity total_instructions = {"TOTAL_INSTRUCTIONS",PAPI_TOT_INS, "PAPI_TOT_INS"};
	static hardwareQuantity flpt_instructions = {"FLPT_INSTRUCTIONS",PAPI_FP_INS, "PAPI_FP_INS"};
    static hardwareQuantity fp_ops = {"FP_OPS", PAPI_FP_OPS, "PAPI_FP_OPS"};

	//Set up the map of all hardware quantities
	 AllHardwareQuantitiesMap = {
			{L1_CACHE_MISS,l1_cache_miss},
			{L2_CACHE_MISS,l2_cache_miss},
			{L3_CACHE_MISS,l3_cache_miss},
			{BRANCH_MISPRED,branch_mispredictions},
			{TOTAL_CYCLES,total_cycles},
			{TOTAL_INSTRUC,total_instructions},
			{FLPT_INSTRUC,flpt_instructions},
            {FP_OPS,fp_ops}
	};

	for (unsigned i = 0; i < counterQuantities.size(); i++)
	{
		std::unordered_map<int,hardwareQuantity>::const_iterator search =
				AllHardwareQuantitiesMap.find(counterQuantities.at(i));
		hardwareQuantitiesMap.insert( {search->first, search->second} );
	}

	for (unsigned i = 0; i < quantities.size(); i++)
	{
		std::unordered_map<int,hardwareQuantity>::const_iterator search =
				AllHardwareQuantitiesMap.find(quantities.at(i));
		hardwareQuantitiesMap.insert( {search->first, search->second} );
	}


}

std::vector<double> GPTLHardwareCounter::getValues() const {

	std::vector<double> papiValues;

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
	for (auto it = hardwareQuantitiesMap.begin(); it != hardwareQuantitiesMap.end(); ++it)
	{
		double papiVal = 0.0;
		std::cout << "papistring = " << it->second.papiQuantityString.c_str() << std::endl;
		GPTLget_eventvalue(this->getName().c_str(), (it->second.papiQuantityString).c_str(), -1, &papiVal);
		papiValues.push_back (papiVal);
	}

	return papiValues;
}


std::vector<std::string> GPTLHardwareCounter::getHardwareQuantities() const {

	std::vector<std::string> namesOfQuantities;

	for (auto it = hardwareQuantitiesMap.begin(); it != hardwareQuantitiesMap.end(); ++it)
	{
		namesOfQuantities.push_back (it->second.hardwareQuantityString);
	}

	return namesOfQuantities;
}



