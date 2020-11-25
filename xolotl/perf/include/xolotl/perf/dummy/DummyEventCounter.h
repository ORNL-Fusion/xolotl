#ifndef DUMMYEVENTCOUNTER_H
#define DUMMYEVENTCOUNTER_H

#include <string>

#include <xolotl/perf/IEventCounter.h>

namespace xolotl
{
namespace perf
{
namespace dummy
{
/**
 * The DummyEventCounter class is instantiated by the DummyHandlerRegistry
 * class and realizes the DummyEventCounter interface.
 */
class DummyEventCounter : public IEventCounter
{
private:
	/**
	 * The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
	DummyEventCounter(void)
	{
	}

public:
	/**
	 * DummyEventCounter constructor that takes the argument name but
	 * doesn't do anything with it
	 */
	DummyEventCounter(const std::string& name)
	{
	}

	/**
	 * The destructor
	 */
	virtual ~DummyEventCounter()
	{
	}

	/**
	 * This operation returns the value of the DummyEventCounter,
	 * the frequency of the specified event.
	 */
	virtual IEventCounter::ValType
	getValue() const
	{
		return 0;
	}

	/**
	 * This operation increments the DummyEventCounter.
	 */
	virtual void
	increment()
	{
	}
};
// end class DummyEventCounter

} // end namespace dummy
} // end namespace perf
} // end namespace xolotl

#endif
