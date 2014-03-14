#ifndef HWQUANTITYINFOMAP_H
#define HWQUANTITYINFOMAP_H

#include <string>
#include <map>
#include "HardwareQuantities.h"

namespace xolotlPerf {

struct HardwareQuantityInfo
{
    // our name for the hardware counter
    std::string name;   

    // PAPI ID code for the hardware counter
    int papiID;

    // PAPI name for the hardware counter
    std::string papiName;

    HardwareQuantityInfo( std::string _name,
                            int _papiID,
                            std::string _papiName )
      : name( _name ),
        papiID( _papiID ),
        papiName( _papiName )
    { }
};


class HardwareQuantityInfoMap : public std::map<HardwareQuantities,HardwareQuantityInfo>
{
public:
    // Initialize the map
    HardwareQuantityInfoMap( void );
};

} // end namespace xolotlPerf

#endif // HWQUANTITYINFOMAP_H

