#include "DummyHandlerRegistry.h"

namespace xolotlPerf
{

// Obtain a Timer by name.
std::shared_ptr<ITimer>
DummyHandlerRegistry::getTimer( std::string name )
{
    // TODO is there a need for us to retain access to this Timer?
    // TODO do we need to check whether client has already created
    // an object with this name and return that object?
    return std::make_shared<DummyTimer>( name );
}


// Obtain an EventCounter by name.
std::shared_ptr<IEventCounter>
DummyHandlerRegistry::getEventCounter( std::string name )
{
    // TODO is there a need for us to retain access to this Timer?
    // TODO do we need to check whether client has already created
    // an object with this name and return that object?
    return std::make_shared<DummyEventCounter>( name );
}

// Obtain a HardwareCounter object by name and by the
// counter data it collects.
std::shared_ptr<IHardwareCounter>
DummyHandlerRegistry::getHardwareCounter( std::string name, 
                                        std::vector<HardwareQuantities> hwq )
{
    // TODO is there a need for us to retain access to this Timer?
    // TODO do we need to check whether client has already created
    // an object with this name and return that object?
    return std::make_shared<DummyHardwareCounter>( name, hwq );
}

// Output any collected performance data to the given output stream.
void DummyHandlerRegistry::dump( std::ostream& /* os */ ) const
{
    // do nothing
}

void DummyHandlerRegistry::dump( int rank) const
{
    // do nothing
}

};  // end namespace xolotlPerf

