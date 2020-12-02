#pragma once

#include <map>
#include <string>

#include <xolotl/perf/IHardwareCounter.h>

namespace xolotl
{
namespace perf
{
namespace papi
{
//! A collection of hardware performance counters.
class PAPIHardwareCounter : public IHardwareCounter
{
private:
	//! Information about a single hardware counter.
	struct CounterSpecInfo
	{
		std::string name; ///< Common name for the counter.
		std::string papiName; ///< PAPI's name for the counter.
		int papiEventID; ///< PAPI event code for the counter.

		CounterSpecInfo(const std::string& _name, const std::string& _papiName,
			int _papiEventID) :
			name(_name),
			papiName(_papiName),
			papiEventID(_papiEventID)
		{
		}
	};

	//! A nicer name for a collection of CounterSpecInfo structs,
	//! keyed by the counter spec ID.
	typedef std::map<IHardwareCounter::CounterSpec, CounterSpecInfo>
		CounterSpecMap;

	//! A collection of CounterSpecInfo structs indicating the
	//! set of hardware counters that we know how to monitor.
	static CounterSpecMap csMap;

	//! The hardware performance counter values we have collected.
	//! These are only valid after the collection has stopped counting.
	IHardwareCounter::ValType vals;

	//! Our configuration (which hardware performance counters
	//! we are monitoring).
	IHardwareCounter::SpecType spec;

	//! Our PAPI event set handle.
	int eventSet;

	//! Construct a PAPIHardwareCounter.
	//! The default constructor is private to force callers to
	//! provide a collection of hardware counters to monitor.
	PAPIHardwareCounter()
	{
	}

	//! Initialize our collection of specifications for the
	//! hardware counters we know how to monitor.
	static void
	initCounterSpecMap();

public:
	/**
	 * Construct a PAPIHardwareCounter.
	 *
	 * @param cset The collection of hardware counter spec values indicating
	 * the set of hardware counters we should monitor
	 */
	PAPIHardwareCounter(const IHardwareCounter::SpecType& cset);

	/// Destroy the counter set.
	virtual ~PAPIHardwareCounter();

	/**
	 * \see IHardwareCounter.h
	 */
	void
	start() override;

	/**
	 * \see IHardwareCounter.h
	 */
	void
	stop() override;

	/**
	 * \see IHardwareCounter.h
	 */
	const ValType&
	getValues() const override
	{
		return vals;
	}

	/**
	 * \see IHardwareCounter.h
	 */
	const SpecType&
	getSpecification() const override
	{
		return spec;
	}

	/**
	 * \see IHardwareCounter.h
	 */
	std::string
	getCounterName(IHardwareCounter::CounterSpec cs) const override;

	/**
	 * Add the given HardwareCounter's value to my value.
	 *
	 * @param cset The counter set whose values should be added to my values.
	 * This counter set must be configured exactly the same as me,
	 * or the operation will throw a std::invalid_argument exception.
	 * @return Myself after adding the given counter set's values.
	 */
	virtual IHardwareCounter&
	operator+=(const IHardwareCounter& c);
};
} // namespace papi
} // namespace perf
} // namespace xolotl
