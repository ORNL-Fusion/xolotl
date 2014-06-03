#include "XolotlConfigPerf.h"
#include "xolotlPerf.h"
#include <DummyHandlerRegistry.h>
#include <iostream>

#if defined(HAVE_PERFLIB_STD)
#include <StandardHandlerRegistry.h>
#endif // defined(HAVE_PERFLIB_STD)

namespace xolotlPerf
{

static std::shared_ptr<IHandlerRegistry> theHandlerRegistry;


// Create the desired type of handler registry.
bool
initialize( bool useStdRegistry,
                std::vector<HardwareQuantities> hwQuantities )
{
    bool ret = true;

    if( useStdRegistry )
    {
#if defined(HAVE_PERFLIB_STD)
        // we are to use a standard handler registry
        // (one that collects timings)
        theHandlerRegistry = std::make_shared<StandardHandlerRegistry>( hwQuantities );
#else
        // TODO is there another mechanism for writing errors
        // e.g., one that logs error messages?
        throw std::string("\nxolotlPerf::initialize: unable to build requested standard handler registry due to missing dependencies");
#endif // defined(HAVE_PERFLIB_STD)
    }
    else
    {
        // use a dummy HandlerRegistry for this run
        // Note that the dummy (stub) handlers don't take the 
        // collection of hardware quantities to monitor, since
        // they don't monitor anything.
        theHandlerRegistry = std::make_shared<xolotlPerf::DummyHandlerRegistry>();
    }

    return ret;
}

// Provide access to our handler registry.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry( void )
{
    if( !theHandlerRegistry )
    {
        // Throw an error since we have not yet been initialized
        throw std::string("\nxolotlPerf handler registry requested, but library has not been initialized");
    }
    return theHandlerRegistry;
}


} // end namespace xolotlPerf

