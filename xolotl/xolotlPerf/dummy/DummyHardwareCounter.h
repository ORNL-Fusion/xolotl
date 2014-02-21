#ifndef DUMMYHARDWARECOUNTER_H
#define DUMMYHARDWARECOUNTER_H

#include <string>
#include <vector>
#include <memory>
#include "HardwareQuantities.h"
#include "HardwareCounter.h"

namespace xolotlPerf{

/**
 * The DummyHardwareCounter class is instantiated by the DummyHandlerRegistry class
 * and realizes the DummyHardwareCounter interface.
 */
class DummyHardwareCounter : public HardwareCounter
{

private:

	/**
	 * The default constructor is private because HardwareCounters must
	 * always be initialized with a name.
	 */
//	DummyHardwareCounter():HardwareCounter() {}

public:

	/**
	 * DummyHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param aname The DummyHardwareCounter's name
	 * @param hquantities The DummyHardwareCounter's list of quantities
	 */
	DummyHardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities);

	/**
	 * The copy constructor.
	 * @param other The DummyHardwareCounter to copy
	 */
//	DummyHardwareCounter(const DummyHardwareCounter &other);

	/**
	 * The destructor
	 */
	~DummyHardwareCounter();

	/**
	 * This operation returns a DummyHardwareCounter that is created using the copy
	 * constructor. If this DummyHardwareCounter is actually a subclass of DummyHardwareCounter, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * DummyHardwareCounter.
	 * @return A copy of this DummyHardwareCounter.
	 */
//	virtual std::shared_ptr<HardwareCounter> clone();

    /**
     * This operation returns a list of values of the, initially specified,
     * PAPI preset quantities monitored by the DummyHardwareCounter.
     */
    std::vector<long long> getValues() const;

    // This operation returns the name of the DummyHardwareCounter.
    const std::string getName() const;

    // This operation increments the DummyHardwareCounter.
    void increment();

};  //end class DummyHardwareCounter

}  //end namespace xolotlPerf

#endif
