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
/// A collection of hardware performance counters.
class PAPIHardwareCounter : public IHardwareCounter
{
private:
	/// Information about a single hardware counter.
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

	/// A nicer name for a collection of CounterSpecInfo structs,
	/// keyed by the counter spec ID.
	typedef std::map<IHardwareCounter::CounterSpec, CounterSpecInfo*>
		CounterSpecMap;

	/// A collection of CounterSpecInfo structs indicating the
	/// set of hardware counters that we know how to monitor.
	static CounterSpecMap csMap;

	/// The hardware performance counter values we have collected.
	/// These are only valid after the collection has stopped counting.
	IHardwareCounter::ValType vals;

	/// Our configuration (which hardware performance counters
	/// we are monitoring).
	IHardwareCounter::SpecType spec;

	/// Our PAPI event set handle.
	int eventSet;

	/// Construct a PAPIHardwareCounter.
	/// The default constructor is private to force callers to
	/// provide a name and a collection of hardware counters to monitor.
	PAPIHardwareCounter()
	{
	}

	/// Initialize our collection of specifications for the
	/// hardware counters we know how to monitor.
	static void
	initCounterSpecMap();

public:
	/// Construct a PAPIHardwareCounter.
	///
	/// @param name The name to associate with the collected counts.
	/// @param cset The collection of hardware counter spec values indicating
	///             The set of hardware counters we should monitor.
	PAPIHardwareCounter(
		const std::string& name, const IHardwareCounter::SpecType& cset);

	/// Destroy the counter set.
	virtual ~PAPIHardwareCounter();

	/// Start counting hardware counter events.
	void
	start() override;

	/// Stop counting hardware counter events.
	void
	stop() override;

	///
	/// Retrieve the values of the hardware counters that have been collected.
	/// The values are only valid if the counter set is not currently counting.
	///
	/// @return The current counts for our configured values.
	///
	const ValType&
	getValues() const override
	{
		return vals;
	}

	///
	/// Retrieve the configuration of the IHardwareCounter.
	///
	/// @return The hardware counters the counter set was configured to collect.
	///
	const SpecType&
	getSpecification() const override
	{
		return spec;
	}

	/// Retrieve the name of the given hardware counter.
	/// @return The name of the given hardware counter.
	std::string
	getCounterName(IHardwareCounter::CounterSpec cs) const override;

	/// Add the given HardwareCounter's value to my value.
	/// @param cset The counter set whose values should be added to my values.
	///             This counter set must be configured exactly the same as me,
	///             or the operation will throw a std::invalid_argument
	///             exception.
	/// @return Myself after adding the given counter set's values.
	virtual IHardwareCounter&
	operator+=(const IHardwareCounter& c);
};
} // namespace papi
} // namespace perf
} // namespace xolotl
