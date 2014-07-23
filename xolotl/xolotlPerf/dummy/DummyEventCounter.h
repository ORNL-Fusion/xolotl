#ifndef DUMMYEVENTCOUNTER_H
#define DUMMYEVENTCOUNTER_H

#include <string>
#include "xolotlCore/Identifiable.h"
#include "xolotlPerf/IEventCounter.h"

using namespace std;

namespace xolotlPerf{

/**
 * The DummyEventCounter class is instantiated by the DummyHandlerRegistry
 * class and realizes the DummyEventCounter interface.
 */
class DummyEventCounter : public IEventCounter, public xolotlCore::Identifiable
{

private:

	/**
	 * The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
    DummyEventCounter(void)
      : xolotlCore::Identifiable("unused")
    { }


public:

	/**
	 * DummyEventCounter constructor that takes the argument name
	 *
	 * @param name The DummyEventCounter's name
	 */
	DummyEventCounter(std::string name)
      : xolotlCore::Identifiable("unused")
    { }


	/**
	 * The destructor
	 */
	virtual ~DummyEventCounter() { }

	/**
	 * This operation returns the value of the DummyEventCounter,
	 * the frequency of the specified event.
	 */
	virtual IEventCounter::ValType getValue() const  { return 0; }


	/**
	 * This operation increments the DummyEventCounter.
	 */
	virtual void increment()    {  }


};  //end class DummyEventCounter

}  //end namespace xolotlPerf

#endif
