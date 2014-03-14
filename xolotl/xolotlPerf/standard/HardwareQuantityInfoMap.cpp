#include "papi.h"
#include "HardwareQuantityInfoMap.h"


namespace xolotlPerf {

HardwareQuantityInfoMap::HardwareQuantityInfoMap( void )
{
    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( L1_CACHE_MISS, HardwareQuantityInfo( "L1_CACHE_MISS", PAPI_L1_TCM, "PAPI_L1_TCM" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( L2_CACHE_MISS, HardwareQuantityInfo( "L2_CACHE_MISS", PAPI_L2_TCM, "PAPI_L2_TCM" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( L3_CACHE_MISS, HardwareQuantityInfo( "L3_CACHE_MISS", PAPI_L3_TCM, "PAPI_L3_TCM" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( BRANCH_MISPRED, HardwareQuantityInfo( "BRANCH_MISPREDICTIONS", PAPI_BR_MSP, "PAPI_BR_MSP" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( TOTAL_CYCLES, HardwareQuantityInfo( "TOTAL_CYCLES", PAPI_TOT_CYC, "PAPI_TOT_CYC" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( TOTAL_INSTRUC, HardwareQuantityInfo( "TOTAL_INSTRUCTIONS", PAPI_TOT_INS, "PAPI_TOT_INS" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( FLPT_INSTRUC, HardwareQuantityInfo( "FLPT_INSTRUCTIONS", PAPI_FP_INS, "PAPI_FP_INS" ) ) );

    insert( std::pair<HardwareQuantities,HardwareQuantityInfo>( FP_OPS, HardwareQuantityInfo( "FP_OPS", PAPI_FP_OPS, "PAPI_FP_OPS" ) ) );

}

} // end namespace xolotlPerf

