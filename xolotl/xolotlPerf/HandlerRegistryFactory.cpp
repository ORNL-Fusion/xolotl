#include "XolotlConfigPerf.h"
#include "HandlerRegistryFactory.h"
#include "dummy/DummyHandlerRegistry.h"
#include <iostream>

#if defined(HAVE_PERFLIB_STD)
#include "standard/StandardHandlerRegistry.h"
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
        std::cerr << "xolotlPerf::initialize: unable to build requested standard handler registry" << std::endl;
        ret = false;
#endif // defined(HAVE_PERFLIB_STD)
    }
    else
    {
        // use a dummy HandlerRegistry for this run
        // Note that the dummy (stub) handlers don't take the 
        // collection of hardware quantities to monitor, since
        // they don't monitor anything.
        theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
    }

    return ret;
}




// Provide access to our handler registry.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry( void )
{
    if( !theHandlerRegistry )
    {
        // We have not yet been initialized.
        // Issue a warning and use a dummy (stub) registry.
        std::cerr << "Warning: xolotlPerf handler registry requested, but library has not been initialized; using dummy handlers" << std::endl;

        xolotlPerf::initialize( false );
    }
    return theHandlerRegistry;
}


} // end namespace xolotlPerf

