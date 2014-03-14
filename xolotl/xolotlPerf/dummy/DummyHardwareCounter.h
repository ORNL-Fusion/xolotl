#ifndef DUMMYHARDWARECOUNTER_H
#define DUMMYHARDWARECOUNTER_H

#include <string>
#include <vector>
#include <memory>
#include "Identifiable.h"
#include "IHardwareCounter.h"
#include "HardwareQuantities.h"


namespace xolotlPerf{

/**
 * The DummyHardwareCounter class is instantiated by the DummyHandlerRegistry class
 * and realizes the DummyHardwareCounter interface.
 */
class DummyHardwareCounter : public IHardwareCounter, public xolotlCore::Identifiable
{

private:

	/**
	 * The default constructor is private because HardwareCounters
	 * must always be given a name and a vector of quantities to
	 * be monitored.
	 */
	DummyHardwareCounter()
		: xolotlCore::Identifiable("unused")
	{ }

public:

	/**
	 * DummyHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param counterName The DummyHardwareCounter's name
	 * @param counterQuantities The vector of quantities the DummyHardwareCounter will monitor
	 */
	DummyHardwareCounter(std::string counterName,
			const std::vector<HardwareQuantities> &counterQuantities)
		: xolotlCore::Identifiable("unused")
	{ }

	/**
	 * The destructor
	 */
	virtual ~DummyHardwareCounter() { }

    /**
     * This operation returns a list of values of the, initially specified,
     * hardware quantities monitored by the DummyHardwareCounter.
     */
    std::vector<double> getValues(void) const;


	/**
	 * This operation returns the list of hardware
	 * quantities monitored by the GPTLHardwareCounter.
	 */
	virtual std::vector<std::string> getHardwareQuantities() const;

};  //end class DummyHardwareCounter

}  //end namespace xolotlPerf

#endif
