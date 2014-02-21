#ifndef DUMMYEVENTCOUNTER_H
#define DUMMYEVENTCOUNTER_H

#include <string>
#include <memory>
#include "EventCounter.h"

using namespace std;

namespace xolotlPerf{

// The DummyEventCounter class is instantiated by the DummyHandlerRegistry
// class and realizes the DummyEventCounter interface.
class DummyEventCounter : public EventCounter
{

private:

	/** The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
	DummyEventCounter():EventCounter("") {}


public:

	/**
	 * DummyEventCounter constructor that takes the argument name
	 *
	 * @param aname The DummyEventCounter's name
	 */
	DummyEventCounter(std::string aname);

	/**
	 * The copy constructor.
	 * @param other The DummyEventCounter to copy
	 */
	DummyEventCounter(const DummyEventCounter &other);

	/** The destructor
	 */
	~DummyEventCounter();

	/**
	 * This operation returns a DummyEventCounter that is created using the copy
	 * constructor. If this DummyEventCounter is actually a subclass of DummyEventCounter, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * DummyEventCounter.
	 * @return A copy of this DummyEventCounter.
	 */
//	virtual std::shared_ptr<EventCounter> clone();

	/**
	 * This operation returns the value of the DummyEventCounter,
	 * the frequency of the specified event.
	 */
	int getValue();

	/**
	 * This operation returns the name of the DummyEventCounter.
	 *
	 * @return The name of this DummyEventCounter
	 */
	const std::string getName() const;

	/**
	 * This operation sets the name of the DummyEventCounter to the
	 * specified name.
	 *
	 * @param name The new name
	 */
//	void setName(std::string aname);

	//This operation increments the DummyEventCounter.
	void increment();


};  //end class DummyEventCounter

}  //end namespace xolotlPerf

#endif
