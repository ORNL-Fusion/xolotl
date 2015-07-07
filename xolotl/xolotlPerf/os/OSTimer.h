#ifndef OSTIMER_H
#define OSTIMER_H

#include "xolotlPerf/perfConfig.h"
#include "xolotlPerf/ITimer.h"
#include "xolotlCore/Identifiable.h"

#if defined(HAVE_CLOCK_GETTIME)
#  include "xolotlPerf/os/CGTTimestamp.h"
#elif defined(HAVE_GETTIMEOFDAY)
#  include "xolotlPerf/os/GTODTimestamp.h"
#else
#  error "No supported timing source is available."
#endif // defined(HAVE_CLOCK_GETTIME)


namespace xolotlPerf {

/// A timer that measures how long something takes to execute.
/// Uses an operating system/runtime timer interface.
class OSTimer : public ITimer, public xolotlCore::Identifiable
{
private:
    /// The timer's value.
    ITimer::ValType val;

    /// When the timer was started.
    /// Will be OSRTimestamp::invalidValue if timer is not running.
    OSRTimestamp startTime;


    /// Construct a timer.
    /// The default constructor is private to force callers to provide a name for the timer object.
    OSTimer( void )
      : xolotlCore::Identifiable( "unused" ),
        val( 0 )
    { }
public:
    ///
    /// Construct a timer.
    ///
    /// @param name The name to associate with the timer.
    OSTimer( const std::string& name )
      : xolotlCore::Identifiable( name ),
        val( 0 ),
        startTime( OSRTimestamp::invalidValue )
    { }


    ///
    /// Destroy the timer.
    ///
    virtual ~OSTimer( void ) { }


    ///
    /// Start the timer.
    /// Throws std::runtime_error if starting a timer that was already started.
    ///
    virtual void start( void );


    ///
    /// Stop the timer.
    /// Throws std::runtime_error if stopping a timer that was not running.
    ///
    virtual void stop( void );


    ///
    /// Reset the timer's value.
    /// Throws std::runtime_error if resetting a timer that was running.
    ///
    virtual void reset( void );



    ///
    /// Determine if the Timer is currently running.
    /// 
    /// @return true if the Timer is running, false otherwise.
    ///
    virtual bool isRunning( void ) const    { return (startTime != OSRTimestamp::invalidValue); }
    ///
    /// Retrieve the value of the timer.
    /// The value is only valid if the timer is not running.
    ///
    /// @return The elapsed time measured by this timer.
    ///
    virtual ITimer::ValType  getValue( void ) const { return val; }


    ///
    /// Retrieve the Timer value's units.
    /// @return The units in which the timer's value is given.
    ///
    virtual std::string getUnits( void ) const;


    /// Add the given Timer's value to my value.
    /// @param t The timer whose value should be added to my value.
    /// @return Myself after adding the given timer's value.
    virtual ITimer& operator+=( const ITimer& t )
    {
        val += t.getValue();
        return *this;
    }
};

} // namespace xolotlPerf

#endif // OSTIMER_H
