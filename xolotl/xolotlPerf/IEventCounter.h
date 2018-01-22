#ifndef IEVENTCOUNTER_H
#define IEVENTCOUNTER_H

#include "mpi.h"
#include <limits.h>
#include "../xolotlCore/IIdentifiable.h"

namespace xolotlPerf {

/**
 * Realizations of this interface are responsible for the collection
 * of event performance counter data.
 */
class IEventCounter: public virtual xolotlCore::IIdentifiable {

public:
	/**
	 * Type of the counter.
	 */
	typedef unsigned long ValType;

	/**
	 * MPI type used when transmitting a ValType.
	 */
	static constexpr MPI_Datatype MPIValType = MPI_UNSIGNED_LONG;

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

}//end namespace xolotlPerf

#endif
