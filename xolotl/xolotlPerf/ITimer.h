#ifndef ITIMER_H
#define ITIMER_H

#include <string>
#include "IIdentifiable.h"

using namespace std;

namespace xolotlPerf {

/**
 * Realizations of this interface are responsible for the collection
 * of performance timing statistics.
 */
class ITimer : public virtual xolotlCore::IIdentifiable {

public:

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
    virtual double getValue() const = 0;

    /**
     * This operation returns the units of the ITimer.
     */
    virtual std::string getUnits() const = 0;

};
//end class ITimer

}//end namespace xolotlPerf

#endif
