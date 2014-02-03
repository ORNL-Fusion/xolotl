#ifndef HARDWAREQUANTITIES_H
#define HARDWAREQUANTITIES_H
//Begin section for file hardwarequantities.h
//TODO: Add definitions that you want preserved
//End section for file hardwarequantities.h


//<p>Enumeration to represent the hardware quantities, to be monitored by HardwareCounter, which will be mapped to their corresponding PAPI counterparts.  </p><p>The PAPI preset events which will be included are PAPI_L1_TCM, PAPI_L2_TCM, PAPI_L3_TCM for the levels 1-3 total cache misses, PAPI_BR_MSP for the number mispredicted branches, PAPI_TOT_CYC and PAPI_TOT_INS for the total number of instructions and cycles executed which can be used for the calculation of the total number of instructions (cycles) per cycle (instruction), and PAPI_FP_INS for the number of floating point instructions executed which can be used for the calculation of floating point operations.</p>

//@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
enum HardwareQuantities
{



        //<p>The number of Level 1 cache misses.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        L1_CACHE_MISS,


        //<p>The number of Level 2 cache misses.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        L2_CACHE_MISS,


        //<p>The number of Level 3 cache misses.</p><p>Misses at this level result in data retrieval from main memory.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        L3_CACHE_MISS,


        //<p>The number of branch mispredictions.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        BRANCH_MISPRED,


        //<p>The total number of cycles executed.</p><p>This quantity can be used with TOTAL_INSTRUC to determine the total number of cycles (instructions) per instruction (cycle).</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TOTAL_CYCLES,


        //<p >The total number of instructions executed.</p><p >This quantity can be used with TOTAL_CYCLES to determine the total number of instructions (cycles) per cycle (instruction).</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        TOTAL_INSTRUC,


        //<p>The number of floating point instructions executed.  </p><p>This quantity is used to compute the number of floating point operations.</p>
        //@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"
        FLPT_INSTRUC

}; //end enum HardwareQuantities 


#endif
