#ifndef DUMMYHARDWARECOUNTER_H
#define DUMMYHARDWARECOUNTER_H

#include <memory>
#include <string>
#include <vector>

#include <xolotl/perf/IHardwareCounter.h>

namespace xolotl
{
namespace perf
{
namespace dummy
{
/**
 * The DummyHardwareCounter class is instantiated by the DummyHandlerRegistry
 * class and realizes the DummyHardwareCounter interface.
 */
class DummyHardwareCounter : public IHardwareCounter
{
private:
	/**
	 * An empty collection of values, used by getValues to avoid
	 * returning a reference to a temporary.
	 */
	static IHardwareCounter::ValType dummyVals;

	/**
	 * An empty HardwareCounter specification, used by getSpecification
	 * to avoid returning a reference to a temporary.
	 */
	static IHardwareCounter::SpecType dummySpec;

	/**
	 * The default constructor is private because HardwareCounters
	 * must always be given a name and a vector of quantities to
	 * be monitored.
	 */
	DummyHardwareCounter()
	{
	}

public:
	/**
	 * DummyHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param counterQuantities The vector of quantities the
	 * DummyHardwareCounter will monitor
	 */
	DummyHardwareCounter(const IHardwareCounter::SpecType& ctrSpec)
	{
	}

	/**
	 * The destructor
	 */
	virtual ~DummyHardwareCounter()
	{
	}

	/**
	 * \see IHardwareCounter.h
	 */
	virtual void
	start(void)
	{
		// nothing to do - this is a stub
	}

	/**
	 * \see IHardwareCounter.h
	 */
	virtual void
	stop(void)
	{
		// nothing to do - this is a stub
	}

	/**
	 * \see IHardwareCounter.h
	 */
	virtual const ValType&
	getValues(void) const
	{
		// nothing to do - this is a stub
		return dummyVals;
	}

	/**
	 * \see IHardwareCounter.h
	 */
	virtual const SpecType&
	getSpecification(void) const
	{
		return dummySpec;
	}

	/**
	 * \see IHardwareCounter.h
	 */
	virtual std::string
	getCounterName(IHardwareCounter::CounterSpec cs) const
	{
		return std::string();
	}
};

} // end namespace dummy
} // end namespace perf
} // end namespace xolotl

#endif
