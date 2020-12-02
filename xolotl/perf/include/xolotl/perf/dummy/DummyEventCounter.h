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
public:
	/**
	 * Default constructor
	 */
	DummyEventCounter()
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
