#ifndef TIMESTAMP_H
#define TIMESTAMP_H

#include <time.h>

namespace xolotlPerf
{

struct OSRTimestamp : public timespec
{
public:
    /// An invalid value for a OSRTimestamp.
    static OSRTimestamp invalidValue;


    /// Construct an OSRTimestamp to an invalid value.
    OSRTimestamp( void )
    {
        this->tv_sec = -1;
        this->tv_nsec = -1;
    }



    /// Convert the OSRTimestamp to a quantity in seconds.
    /// @return The timestamp's value as a quantity in seconds.
    ITimer::ValType ToSeconds( void ) const
    {
        return this->tv_sec + (this->tv_nsec / 1.0e9);
    }


    /// Set the OSRTimestamp to the current time.
    void SetToCurrentTime( void )
    {
        clock_gettime( CLOCK_REALTIME, this );
    }


    /// Determine if two OSRTimestamps are equal.
    /// @param ts A timestamp to test for equality.
    /// @return true if the given timestamp is equal to this one, false otherwise.
    bool operator==( const OSRTimestamp& ts ) const
    {
        return (this->tv_sec == ts.tv_sec) && (this->tv_nsec == ts.tv_nsec);
    }

    
    /// Determine if two OSRTimestamps are not equal.
    /// @param ts A timestamp to test for inequality.
    /// @return true if the given timestamp is not equal to this one, false otherwise.
    bool operator!=( const OSRTimestamp& ts ) const
    {
        return !this->operator==(ts);
    }


    /// Subtract the given OSRTimestamp from our value and update our value
    /// with the difference.
    /// @param ts The timestamp to subtract.
    /// @return Ourself after subtracting the given timestamp.
    OSRTimestamp& operator-=( const OSRTimestamp& ts )
    {
        if( ts.tv_nsec > this->tv_nsec )
        {
            // We would have underflow if we subtracted the usec members as is.
            // Adjust the values so that we won't underflow when we subtract.
            this->tv_sec -= 1;
            this->tv_nsec += 1e9;
        }
        this->tv_sec -= ts.tv_sec;
        this->tv_nsec -= ts.tv_nsec;

        return *this;
    }
};

/// Subtract OSRTimestamp b from a.
/// @param a The minuend.
/// @param b The subtrahend.
/// @return The difference between a and b (a - b).
inline
OSRTimestamp 
operator-(OSRTimestamp a, const OSRTimestamp& b)
{
    // Note: parameter a is passed by value, parameter b by reference, 
    // so caller's version of a will not be changed.
    return a -= b;
}

} // namespace xolotlPerf

#endif // TIMESTAMP_H
