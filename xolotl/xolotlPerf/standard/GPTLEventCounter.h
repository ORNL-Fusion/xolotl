#ifndef GPTLEVENTCOUNTER_H
#define GPTLEVENTCOUNTER_H
#include "/home/cxj/Libraries/GPTL/GPTL/gptl-5.0/include/gptl.h"
#include "EventCounter.h"

namespace xolotlPerf{

/**
 * The GPTLEventCounter class is instantiated by the StandardHandlerRegistry
 * class and realizes the EventCounter interface to access event performance
 * counter data found via the General Purpose Timing Library (GPTL).
 */
class GPTLEventCounter : public EventCounter
{

private:

	/**
	 * The value of this EventCounter.
	 */
//	long long value;
//	long value;
	int value;

	/**
	 * The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
	GPTLEventCounter():EventCounter("") { }


public:

	/**
	 * GPTLEventCounter constructor that takes the argument name
	 *
	 * @param aname The GPTLEventCounter's name
	 */
	GPTLEventCounter(std::string aname);

	/**
	 * The copy constructor.
	 * @param other The GPTLEventCounter to copy
	 */
	GPTLEventCounter(const GPTLEventCounter &other);

	/**
	 * The destructor
	 */
	~GPTLEventCounter();

	/**
	 * This operation returns a GPTLEventCounter that is created using the copy
	 * constructor. If this GPTLEventCounter is actually a subclass of GPTLEventCounter, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * GPTLEventCounter.
	 * @return A copy of this GPTLEventCounter.
	 */
//	virtual std::shared_ptr<EventCounter> clone();

	/**
	 * This operation returns the value of the GPTLEventCounter,
	 * the frequency of the specified event.
	 */
	int getValue();

	/**
	 * This operation returns the name of the GPTLEventCounter.
	 *
	 * @return The name of this GPTLEventCounter
	 */
	const std::string getName() const;

	/**
	 * This operation sets the name of the GPTLEventCounter to the
	 * specified name.
	 *
	 * @param name The new name
	 */
//	void setName(std::string aname);

	//This operation increments the GPTLEventCounter.
	void increment();

};  //end class GPTLEventCounter


}  //end namespace xolotlPerf

#endif
