#include <float.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <tuple>

#include <xolotl/perf/EventCounter.h>
#include <xolotl/perf/standard/StdHandlerRegistry.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace perf
{
namespace standard
{
StdHandlerRegistry::StdHandlerRegistry(void)
{
	// nothing else to do
}

StdHandlerRegistry::~StdHandlerRegistry(void)
{
	// Release the objects we have been tracking.
	// Because we use shared_ptrs for these objects,
	// we do not need to explicitly delete the objects themselves.
	allTimers.clear();
	allEventCounters.clear();
	allHWCounterSets.clear();
}

// We can create the EventCounters, since they don't depend on
// more specialized functionality from any of our subclasses.
std::shared_ptr<IEventCounter>
StdHandlerRegistry::getEventCounter(const std::string& name)
{
	// TODO - associate the object we create with the current region
	std::shared_ptr<IEventCounter> ret;

	// Check if we have already created an event counter with this name.
	auto iter = allEventCounters.find(name);
	if (iter != allEventCounters.end()) {
		// We have already created an event counter with this name.
		// Return name.
		ret = iter->second;
	}
	else {
		// We have not yet created an event counter with this name.
		// Build one and keep track of it.
		ret = std::make_shared<EventCounter>(name);
		allEventCounters[name] = ret;
	}
	return ret;
}

template <typename T, typename V>
void
StdHandlerRegistry::CollectAllObjectNames(int myRank,
	const std::map<std::string, std::shared_ptr<T>>& myObjs,
	std::map<std::string, PerfObjStatistics<V>>& stats) const
{
	// Collect my own object's names.
	std::vector<std::string> myNames;
	CollectMyObjectNames(myObjs, myNames);
	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();

	// Determine amount of space required for names
	unsigned int nBytes = 0;
	for (auto nameIter = myNames.begin(); nameIter != myNames.end();
		 ++nameIter) {
		// Add enough space for the name plus a NUL terminating character.
		nBytes += (nameIter->length() + 1);
	}

	// Let root know how much space it needs to collect all object names
	unsigned int totalNumBytes = 0;
	MPI_Reduce(
		&nBytes, &totalNumBytes, 1, MPI_UNSIGNED, MPI_SUM, 0, xolotlComm);

	// Marshal all our object names.
	char* myNamesBuf = new char[nBytes];
	char* pName = myNamesBuf;
	for (auto nameIter = myNames.begin(); nameIter != myNames.end();
		 ++nameIter) {
		strcpy(pName, nameIter->c_str());
		pName += (nameIter->length() + 1); // skip the NUL terminator
	}
	assert(pName == (myNamesBuf + nBytes));

	// Provide all names to root.
	// First, provide the amount of data from each process.
	int cwSize;
	MPI_Comm_size(xolotlComm, &cwSize);
	char* allNames = (myRank == 0) ? new char[totalNumBytes] : NULL;
	int* allNameCounts = (myRank == 0) ? new int[cwSize] : NULL;
	int* allNameDispls = (myRank == 0) ? new int[cwSize] : NULL;

	MPI_Gather(&nBytes, 1, MPI_INT, allNameCounts, 1, MPI_INT, 0, xolotlComm);

	// Next, root computes the displacements for data from each process.
	if (myRank == 0) {
		allNameDispls[0] = 0;
		for (int i = 1; i < cwSize; ++i) {
			allNameDispls[i] = allNameDispls[i - 1] + allNameCounts[i - 1];
		}
	}

	// Finally, gather all names to the root process.
	MPI_Gatherv(myNamesBuf, nBytes, MPI_CHAR, allNames, allNameCounts,
		allNameDispls, MPI_CHAR, 0, xolotlComm);

	if (myRank == 0) {
		// Process the gathered names to determine the
		// set of all known object names.
		pName = allNames;
		while (pName < (allNames + totalNumBytes)) {
			auto iter = stats.find(pName);
			if (iter == stats.end()) {
				// This is an object  name we have not seen before.
				// Add it to the statistics map.
				stats.insert(std::pair<std::string, PerfObjStatistics<V>>(
					pName, PerfObjStatistics<V>(pName)));
			}

			// Advance to next object name
			pName += (strlen(pName) + 1);
		}
		assert(pName == allNames + totalNumBytes);
	}

	// clean up
	delete[] myNamesBuf;
	delete[] allNames;
	delete[] allNameCounts;
	delete[] allNameDispls;
}

template <typename T>
void
StdHandlerRegistry::CollectMyObjectNames(
	const std::map<std::string, std::shared_ptr<T>>& myObjs,
	std::vector<std::string>& objNames) const
{
	for (auto oiter = myObjs.begin(); oiter != myObjs.end(); ++oiter) {
		objNames.push_back(oiter->first);
	}
}

// Specialize for hardware counter objects, since a given object name
// represents multiple hardware counters.
template <>
void
StdHandlerRegistry::CollectMyObjectNames<IHardwareCounter>(
	const std::map<std::string, std::shared_ptr<IHardwareCounter>>& myObjs,
	std::vector<std::string>& objNames) const
{
	for (auto oiter = myObjs.begin(); oiter != myObjs.end(); ++oiter) {
		std::string baseName = oiter->first;

		const IHardwareCounter::SpecType& spec =
			oiter->second->getSpecification();
		for (auto siter = spec.begin(); siter != spec.end(); ++siter) {
			std::ostringstream namestr;
			namestr << baseName << ':' << oiter->second->getCounterName(*siter);
			objNames.push_back(namestr.str());
		}
	}
}

template <typename T, typename V>
std::pair<bool, V>
StdHandlerRegistry::GetObjValue(
	const std::map<std::string, std::shared_ptr<T>>& myObjs,
	const std::string& objName) const
{
	auto currObjIter = myObjs.find(objName);
	bool found = currObjIter != myObjs.end();
	V val = V();
	if (found) {
		val = currObjIter->second->getValue();
	}
	else {
		// val's value is undefined.
		// Callers must test the found part of the return pair to
		// know whether the value part is valid.
	}
	return std::make_pair(found, val);
}

template <>
std::pair<bool, IHardwareCounter::CounterType>
StdHandlerRegistry::GetObjValue(
	const std::map<std::string, std::shared_ptr<IHardwareCounter>>& myObjs,
	const std::string& objName) const
{
	// Split the objName into an IHardwareCounter object name and a
	// hardware counter name.
	size_t idx = objName.find_last_of(':');
	assert(idx != std::string::npos);
	std::string objNameOnly = objName.substr(0, idx);
	std::string ctrName = objName.substr(idx + 1);

	auto objIter = myObjs.find(objNameOnly);
	bool found = objIter != myObjs.end();
	IHardwareCounter::CounterType val = 0;
	if (found) {
		const std::shared_ptr<IHardwareCounter>& currObj = objIter->second;

		// We know about the object.
		// Does the object contain data for the requested hardware counter?
		unsigned int idx = 0;
		IHardwareCounter::SpecType spec = currObj->getSpecification();
		for (auto siter = spec.begin(); siter != spec.end(); ++siter, ++idx) {
			if (currObj->getCounterName(*siter) == ctrName) {
				break;
			}
		}

		if (idx != spec.size()) {
			// The current object did collect data for the
			// requested counter.  Retrieve that value.
			val = currObj->getValues()[idx];
		}
	}

	return std::make_pair(found, val);
}

template <typename T, typename V>
void
StdHandlerRegistry::AggregateStatistics(int myRank,
	const std::map<std::string, std::shared_ptr<T>>& myObjs,
	std::map<std::string, PerfObjStatistics<V>>& stats) const
{
	// Determine the set of object names known across all processes.
	// Since some processes may define an object that others don't, we
	// have to form the union across all processes.
	// Unfortunately, because the strings are of different lengths,
	// we have a more difficult marshal/unmarshal problem than we'd like.
	CollectAllObjectNames<T, V>(myRank, myObjs, stats);
	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();

	// Let all processes know how many statistics we will be collecting.
	int nObjs;
	if (myRank == 0) {
		nObjs = stats.size();
	}
	MPI_Bcast(&nObjs, 1, MPI_INT, 0, xolotlComm);
	assert(nObjs >= 0);

	// Collect and compute statistics for each object.
	auto tsiter = stats.begin();
	for (int idx = 0; idx < nObjs; ++idx) {
		// broadcast the current object's name
		int nameLen = (myRank == 0) ? tsiter->second.name.length() : -1;
		MPI_Bcast(&nameLen, 1, MPI_INT, 0, xolotlComm);
		// we can safely cast away const on the tsiter data string because
		// the only process that accesses that string is rank 0,
		// and it only reads the data.
		char* objName = (myRank == 0) ?
			const_cast<char*>(tsiter->second.name.c_str()) :
			new char[nameLen + 1];
		MPI_Bcast(objName, nameLen + 1, MPI_CHAR, 0, xolotlComm);

		// do we know about the current object?
		bool knowObject;
		V myVal;
		std::tie<bool, V>(knowObject, myVal) =
			GetObjValue<T, V>(myObjs, objName);

		// collect count of processes knowing the current object
		unsigned int* pcount =
			(myRank == 0) ? &(tsiter->second.processCount) : NULL;
		int knowObjVal = knowObject ? 1 : 0;
		MPI_Reduce(&knowObjVal, pcount, 1, MPI_INT, MPI_SUM, 0, xolotlComm);

		// collect min value of current object
		V* pMinVal = (myRank == 0) ? &(tsiter->second.min) : NULL;
		V reduceVal = knowObject ? myVal : T::MaxValue;
		MPI_Reduce(
			&reduceVal, pMinVal, 1, T::MPIValType, MPI_MIN, 0, xolotlComm);

		// collect max value of current object
		V* pMaxVal = (myRank == 0) ? &(tsiter->second.max) : NULL;
		reduceVal = knowObject ? myVal : T::MinValue;
		MPI_Reduce(
			&reduceVal, pMaxVal, 1, T::MPIValType, MPI_MAX, 0, xolotlComm);

		// collect sum of current object's values (for computing avg and stdev)
		double valSum;
		// use the same myVal as for max: actual value if known, 0 otherwise
		double myValAsDouble = (double)reduceVal;
		MPI_Reduce(
			&myValAsDouble, &valSum, 1, MPI_DOUBLE, MPI_SUM, 0, xolotlComm);
		if (myRank == 0) {
			tsiter->second.average = valSum / tsiter->second.processCount;
		}

		// collect sum of squares of current object's values (for stdev)
		double valSquaredSum;
		double myValSquared = myValAsDouble * myValAsDouble;
		MPI_Reduce(&myValSquared, &valSquaredSum, 1, MPI_DOUBLE, MPI_SUM, 0,
			xolotlComm);
		if (myRank == 0) {
			tsiter->second.stdev =
				sqrt((valSquaredSum / tsiter->second.processCount) -
					(tsiter->second.average * tsiter->second.average));
		}

		// clean up
		if (myRank != 0) {
			delete[] objName;
		}

		// advance to next object
		if (myRank == 0) {
			++tsiter;
		}
	}
}

void
StdHandlerRegistry::collectStatistics(
	PerfObjStatsMap<ITimer::ValType>& timerStats,
	PerfObjStatsMap<IEventCounter::ValType>& counterStats,
	PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats)
{
	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();
	int myRank;
	MPI_Comm_rank(xolotlComm, &myRank);

	// Aggregate statistics about counters in all processes.
	// First, timers...
	AggregateStatistics<ITimer, ITimer::ValType>(myRank, allTimers, timerStats);

	// ...next event counters...
	AggregateStatistics<IEventCounter, IEventCounter::ValType>(
		myRank, allEventCounters, counterStats);

	// ...finally hardware counters.
	AggregateStatistics<IHardwareCounter, IHardwareCounter::CounterType>(
		myRank, allHWCounterSets, hwCounterStats);
}

void
StdHandlerRegistry::reportStatistics(std::ostream& os,
	const PerfObjStatsMap<ITimer::ValType>& timerStats,
	const PerfObjStatsMap<IEventCounter::ValType>& counterStats,
	const PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats) const
{
	os << "\nTimers:\n";
	for (auto iter = timerStats.begin(); iter != timerStats.end(); ++iter) {
		iter->second.outputTo(os);
	}

	os << "\nCounters:\n";
	for (auto iter = counterStats.begin(); iter != counterStats.end(); ++iter) {
		iter->second.outputTo(os);
	}

	os << "\nHardwareCounters:\n";
	for (auto iter = hwCounterStats.begin(); iter != hwCounterStats.end();
		 ++iter) {
		iter->second.outputTo(os);
	}
}

} // namespace standard
} // namespace perf
} // namespace xolotl
