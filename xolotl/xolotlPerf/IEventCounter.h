#ifndef IEVENTCOUNTER_H
#define IEVENTCOUNTER_H

#include "IIdentifiable.h"

namespace xolotlPerf {

//Realizations of this interface are responsible for the collection
//of event performance counter data.
class IEventCounter : public virtual xolotlCore::IIdentifiable {

//private:

	/**
	 * The default constructor is declared private since all event counters
	 * must be initialized with a name.
	 */
//	IEventCounter();

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

};
//end class IEventCounter

}//end namespace xolotlPerf

#endif
