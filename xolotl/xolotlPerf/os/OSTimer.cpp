#include <stdexcept>
#include "xolotlPerf/os/OSTimer.h"

namespace xolotlPerf
{

OSRTimestamp OSRTimestamp::invalidValue;


void
OSTimer::start(void)
{
    if( isRunning() )
    {
        throw std::runtime_error( "Attempting to start a timer that is already running." );
    }
    
    // Start the timer by sampling the current time.
    startTime.SetToCurrentTime();
}



void
OSTimer::stop(void)
{
    if( !isRunning() )
    {
        throw std::runtime_error( "Attempting to stop a timer that was not running." );
    }

    // Stop the timer by sampling the ending time.
    OSRTimestamp endTimestamp;
    endTimestamp.SetToCurrentTime();

    // Form the difference between the end timestamp and 
    // our saved start timestamp.
    val += (endTimestamp - startTime).ToSeconds();

    // Indicate the timer is no longer running.
    startTime = OSRTimestamp::invalidValue;
}


void
OSTimer::reset(void)
{
    if( isRunning() )
    {
        throw std::runtime_error( "Attempting to reset a timer that is running" );
    }

    val = 0;
}



std::string 
OSTimer::getUnits(void) const
{
    return std::string("s");
}

} // namespace xolotlPerf

