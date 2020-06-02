#ifndef IEVENTCOUNTER_H
#define IEVENTCOUNTER_H

#include <mpi.h>
#include <limits.h>
#include <xolotl/core/IIdentifiable.h>

namespace xolotl {
namespace perf {

/**
 * Realizations of this interface are responsible for the collection
 * of event performance counter data.
 */
class IEventCounter: public virtual core::IIdentifiable {

public:
	/**
	 * Type of the counter.
	 */
	typedef unsigned long ValType;

	/**
	 * MPI type used when transmitting a ValType.
	 */
	static const MPI_Datatype MPIValType;

	/**
	 * The minimum value possible.
	 */
	static const ValType MinValue = 0;

	/**
	 * The maximum value possible.
	 */
	static const ValType MaxValue = ULONG_MAX;

	/**
	 * The destructor
	 */
	virtual ~IEventCounter() {
	}

	/**
	 * This operation returns the value of the IEventCounter, the frequency
	 * of the specified event.
	 */
	virtual ValType getValue() const = 0;

	/**
	 * This operation increments the IEventCounter.
	 */
	virtual void increment() = 0;

};
//end class IEventCounter

}//end namespace perf
}//end namespace xolotl

#endif
