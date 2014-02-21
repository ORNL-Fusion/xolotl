#ifndef EVENTCOUNTER_H
#define EVENTCOUNTER_H

#include <string>
#include <memory>

using namespace std;

namespace xolotlPerf {

//Realizations of this interface are responsible for the collection
//of event performance counter data.
class EventCounter {

private:

	/**
	 * The default constructor is declared private since all timers
	 * must be initialized with a name.
	 */
	EventCounter();

protected:

	/**
	 * The name of this EventCounter.
	 */
	std::string name;

	/**
	 * The value of this EventCounter.
	 */
//	long long value;
//	long value;
//	int value;

public:

	/**
	 * EventCounter constructor that takes the argument name
	 *
	 * @param aname The EventCounter's name
	 */
	EventCounter(std::string aname);

	/**
	 * The copy constructor.
	 * @param other The EventCounter to copy
	 */
	EventCounter(const EventCounter &other);

	/**
	 * The destructor
	 */
	virtual ~EventCounter();

	/**
	 * This operation returns a EventCounter that is created using the copy
	 * constructor. If this EventCounter is actually a subclass of EventCounter, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * EventCounter.
	 * @return A copy of this EventCounter.
	 */
//	virtual std::shared_ptr<EventCounter> clone();

	/**
	 * This operation returns the value of the EventCounter, the frequency
	 * of the specified event.
	 */
	virtual int getValue() = 0;

	/**
	 * This operation returns the name of the EventCounter.
	 *
	 * @return The name of this EventCounter
	 */
	virtual const std::string getName() const;
//	const std::string getName() const {return name;}

	/**
	 * This operation sets the name of the EventCounter to the
	 * specified name.
	 *
	 * @param name The new name
	 */
//	void setName(std::string aname);

	//This operation increments the EventCounter.
	virtual void increment() = 0;

};
//end class EventCounter

}//end namespace xolotlPerf

#endif
