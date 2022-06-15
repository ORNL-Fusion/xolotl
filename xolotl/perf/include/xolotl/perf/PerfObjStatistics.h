#ifndef PERFOBJSTATISTICS_H
#define PERFOBJSTATISTICS_H

#include <iostream>
#include <map>
#include <string>

namespace xolotl
{
namespace perf
{
/**
 * Statistics (min,max,avg,stdev) for a performance metric we collected
 * during a program run.
 */
template <class T>
struct PerfObjStatistics
{
	std::string name; //< name of the metric
	unsigned int processCount; //< number of processes that collected data for
							   // this metric
	T min; //< min value across all processes that collected data for the metric
	T max; //< max value across all processes that collected data for the metric
	double average; //< average value across all processes that collected data
					// for the metric
	double stdev; //< standard deviation of values across all processes that
				  // collected data for the metric

	/**
	 * Construct a PerfObjStatistics struct with default values.
	 * @param _name The metric name.
	 */
	PerfObjStatistics(const std::string& _name) :
		name(_name),
		processCount(0),
		min(0),
		max(0),
		average(0),
		stdev(0)
	{
	}

	/**
	 * Construct a PerfObjStatistics struct as a copy of another.
	 * @param obj The object to be copied.
	 */
	PerfObjStatistics(const PerfObjStatistics& obj) :
		name(obj.name),
		processCount(obj.processCount),
		min(obj.min),
		max(obj.max),
		average(obj.average),
		stdev(obj.stdev)
	{
	}

	/**
	 * Replace my own values to be a copy of another PerfObjStatistics.
	 * @param obj The object to be copied.
	 */
	PerfObjStatistics&
	operator=(const PerfObjStatistics& obj)
	{
		if (&obj != this) {
			name = obj.name;
			processCount = obj.processCount;
			min = obj.min;
			max = obj.max;
			average = obj.average;
			stdev = obj.stdev;
		}
		return *this;
	}

	/**
	 * Output our name and statistics to the given output stream.
	 * @param os The output stream on which we will write our statistics.
	 */
	void
	outputTo(std::ostream& os) const
	{
		// Output data in YAML format.
		const char* nameIndent = "  ";
		const char* propIndent = "    ";
		os << nameIndent << name << ":\n"
		   << propIndent << "process_count: " << processCount << '\n'
		   << propIndent << "min: " << min << '\n'
		   << propIndent << "max: " << max << '\n'
		   << propIndent << "average: " << average << '\n'
		   << propIndent << "stdev: " << stdev << '\n'
		   << std::endl;
	}
};

/* C++11 supports template typedefs.   Earlier C++ standards did not,
 * but we are requiring C++11 for the shared_ptr support so we know
 * that we are using a C++11 compiler.
 */
template <typename T>
using PerfObjStatsMap = std::map<std::string, PerfObjStatistics<T>>;

} // end namespace perf
} // end namespace xolotl

#endif // PERFOBJSTATISTICS_H
