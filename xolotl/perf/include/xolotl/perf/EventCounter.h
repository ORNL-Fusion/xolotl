#ifndef EVENTCOUNTER_H
#define EVENTCOUNTER_H

#include <xolotl/perf/IEventCounter.h>

namespace xolotl
{
namespace perf
{
/**
 * An EventCounter keeps a count.  Code using an EventCounter can
 * increment the counter whenever an event of interest occurs, and
 * retrieve the count whenever it is desired.
 */
class EventCounter : public IEventCounter
{
private:
	/**
	 * The value of this IEventCounter.
	 */
	IEventCounter::ValType value;

public:
	/**
	 * Default constructor
	 */
	EventCounter() : value(0)
	{
	}

	/**
	 * The destructor
	 */
	virtual ~EventCounter()
	{
	}

	/**
	 * This operation returns the value of the EventCounter,
	 * the frequency of the specified event.
	 */
	IEventCounter::ValType
	getValue() const override
	{
		return value;
	}

	/**
	 * This operation increments the EventCounter.
	 */
	void
	increment() override
	{
		++value;
	}

	/**
	 * Add the given counter's value to my value
	 */
	IEventCounter&
	operator+=(const IEventCounter& c) override
	{
		value += c.getValue();
		return *this;
	}

	std::shared_ptr<IEventCounter>
	copy() const override
	{
		return std::make_shared<EventCounter>(*this);
	}
};
// end class EventCounter

} // end namespace perf
} // end namespace xolotl

#endif
