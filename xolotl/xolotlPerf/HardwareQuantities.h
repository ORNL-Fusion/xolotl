#ifndef HARDWAREQUANTITIES_H
#define HARDWAREQUANTITIES_H

namespace xolotlPerf{

// Enumeration to represent the hardware quantities, to be monitored by HardwareCounter, which will be mapped to their corresponding PAPI counterparts.   <p>The PAPI preset events which will be included are PAPI_L1_TCM, PAPI_L2_TCM, PAPI_L3_TCM for the levels 1-3 total cache misses, PAPI_BR_MSP for the number mispredicted branches, PAPI_TOT_CYC and PAPI_TOT_INS for the total number of instructions and cycles executed which can be used for the calculation of the total number of instructions (cycles) per cycle (instruction), and PAPI_FP_INS for the number of floating point instructions executed which can be used for the calculation of floating point operations.

enum HardwareQuantities
{

        // The number of Level 1 cache misses.
        L1_CACHE_MISS,

        // The number of Level 2 cache misses.
        L2_CACHE_MISS,

        // The number of Level 3 cache misses. Misses at this level result in data retrieval from main memory.
        L3_CACHE_MISS,

        // The number of branch mispredictions.
        BRANCH_MISPRED,

        // The total number of cycles executed. This quantity can be used with TOTAL_INSTRUC to determine the total number of cycles (instructions) per instruction (cycle).
        TOTAL_CYCLES,

        //The total number of instructions executed. This quantity can be used with TOTAL_CYCLES to determine the total number of instructions (cycles) per cycle (instruction).
        TOTAL_INSTRUC,

        // The number of floating point instructions executed.  This quantity is used to compute the number of floating point operations.
        FLPT_INSTRUC

}; //end enum HardwareQuantities

}  //end namespace xolotlPerf

#endif
