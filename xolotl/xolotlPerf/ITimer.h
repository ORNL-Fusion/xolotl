#ifndef ITIMER_H
#define ITIMER_H

#include "mpi.h"
#include <string>
#include <float.h>
#include "../xolotlCore/IIdentifiable.h"

using namespace std;

namespace xolotlPerf {

/**
 * Realizations of this interface are responsible for the collection
 * of performance timing statistics.
 */
class ITimer : public virtual xolotlCore::IIdentifiable {

public:
    
    /**
     * The type of a timer value.
     */
    typedef double ValType;


    /**
     * The MPI type to use when transferring a ValType.
     */
    static constexpr MPI_Datatype MPIValType = MPI_DOUBLE;


    /**
     * The minimum value possible.
     */
    static constexpr ValType MinValue = 0.0;


    /**
     * The maximum value possible.
     */
    static constexpr ValType MaxValue = DBL_MAX;


	/**
	 * The destructor
	 */
	virtual ~ITimer() { }


    /**
     * This operations starts the ITimer.
     */
    virtual void start() = 0;

    /**
     * This operation stops the ITimer.
     */
    virtual void stop() = 0;

    /**
     * This operation returns the value of the ITimer.
     */
    virtual ValType getValue() const = 0;

    /**
     * This operation returns the units of the ITimer.
     */
    virtual std::string getUnits() const = 0;

};
//end class ITimer

}//end namespace xolotlPerf

#endif
