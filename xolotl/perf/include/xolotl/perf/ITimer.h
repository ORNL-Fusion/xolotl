#ifndef ITIMER_H
#define ITIMER_H

#include <float.h>

#include <string>

#include <mpi.h>

namespace xolotl
{
namespace perf
{
/**
 * Realizations of this interface are responsible for the collection
 * of performance timing statistics.
 */
class ITimer
{
public:
	/**
	 * The type of a timer value.
	 */
	typedef double ValType;

	/**
	 * The MPI type to use when transferring a ValType.
	 */
	static const MPI_Datatype MPIValType;

	/**
	 * The minimum value possible.
	 */
	static constexpr ValType MinValue = 0.0;

	/**
	 * The maximum value possible.
	 */
	static constexpr ValType MaxValue = DBL_MAX;

	/**
	 * Destroy the timer.
	 */
	virtual ~ITimer()
	{
	}

	/**
	 * Start the timer.
	 */
	virtual void
	start() = 0;

	/**
	 * Stop the timer.
	 */
	virtual void
	stop() = 0;

	/**
	 * Access the timer's value.  (Only valid if timer is not running.)
	 */
	virtual ValType
	getValue() const = 0;

	/**
	 * Reset the timer's value.  Only valid if timer is not running.
	 */
	virtual void
	reset() = 0;

	/**
	 * Obtain a string describing the units of the timer's value.
	 */
	virtual std::string
	getUnits() const = 0;
};
// end class ITimer

} // end namespace perf
} // end namespace xolotl

#endif
