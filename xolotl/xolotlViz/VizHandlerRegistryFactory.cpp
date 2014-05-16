#include "XolotlConfigViz.h"
#include "VizHandlerRegistryFactory.h"
#include <DummyHandlerRegistry.h>
#include <iostream>

#if defined(HAVE_VIZLIB_STD)
#include <StandardHandlerRegistry.h>
#endif // defined(HAVE_VIZLIB_STD)


namespace xolotlViz
{

static std::shared_ptr<IVizHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
bool initialize() {
    bool ret = true;

#if defined(HAVE_VIZLIB_STD)
    // we are to use a standard handler registry
    theHandlerRegistry = std::make_shared<StandardHandlerRegistry>();
#else
    // we are to use a dummy handler registry
    theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
#endif // defined(HAVE_VIZLIB_STD)

    return ret;
}

// Provide access to our handler registry.
std::shared_ptr<IVizHandlerRegistry> getVizHandlerRegistry( void )
{
    if( !theHandlerRegistry )
    {
        xolotlViz::initialize();
    }

    return theHandlerRegistry;
}

} // end namespace xolotlViz

