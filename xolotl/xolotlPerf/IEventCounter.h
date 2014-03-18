#ifndef IEVENTCOUNTER_H
#define IEVENTCOUNTER_H

#include "../xolotlCore/IIdentifiable.h"

namespace xolotlPerf {

/**
 * Realizations of this interface are responsible for the collection
 * of event performance counter data.
 */
class IEventCounter : public virtual xolotlCore::IIdentifiable {

public:

	/**
	 * The destructor
	 */
	virtual ~IEventCounter(){}

	/**
	 * This operation returns the value of the IEventCounter, the frequency
	 * of the specified event.
	 */
	virtual unsigned long getValue() const = 0;

	/**
	 * This operation increments the IEventCounter.
	 */
	virtual void increment() = 0;

}; //end class IEventCounter

} //end namespace xolotlPerf

#endif
