#ifndef DUMMYTIMER_H
#define DUMMYTIMER_H

#include <xolotl/perf/ITimer.h>

namespace xolotl
{
namespace perf
{
namespace dummy
{
/**
 * The DummyTimer class is instantiated by the DummerHandlerRegistry class
 * and realizes the DummyTimer interface.
 */
class DummyTimer : public ITimer
{
public:
	DummyTimer(void)
	{
	}

	virtual ~DummyTimer(void)
	{
	}

	/**
	 * \see ITimer.h
	 */
	virtual void
	start(void)
	{
		return;
	}

	/**
	 * \see ITimer.h
	 */
	virtual void
	stop(void)
	{
		return;
	}

	/**
	 * \see ITimer.h
	 */
	virtual void
	reset(void)
	{
		return;
	}

	/**
	 * \see ITimer.h
	 */
	virtual ITimer::ValType
	getValue(void) const
	{
		return 0;
	}

	/**
	 * \see ITimer.h
	 */
	virtual std::string
	getUnits(void) const
	{
		return "";
	}
};
// end class DummyTimer

} // end namespace dummy
} // end namespace perf
} // end namespace xolotl

#endif
