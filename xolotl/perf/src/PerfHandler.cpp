#include <float.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <sstream>

#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/perf/EventCounter.h>
#include <xolotl/perf/PerfHandler.h>
#include <xolotl/perf/dummy/DummyHardwareCounter.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace perf
{
// OpenMPI 2.x and 3.0 define its MPI_Datatype constants using a
// C-style cast to void*.  Clang++ objects to using these with
// in-class initializers and constexpr in our classes like ITimer.
// So we have to define them elsewhere (i.e., here).
const MPI_Datatype ITimer::MPIValType = MPI_DOUBLE;
const MPI_Datatype IEventCounter::MPIValType = MPI_UNSIGNED_LONG;
const MPI_Datatype IHardwareCounter::MPIValType = MPI_LONG_LONG_INT;

PerfHandler::PerfHandler(const options::IOptions&)
{
}

PerfHandler::~PerfHandler()
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
PerfHandler::getEventCounter(const std::string& name)
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
		ret = std::make_shared<EventCounter>();
		allEventCounters[name] = ret;
	}
	return ret;
}

std::shared_ptr<IHardwareCounter>
PerfHandler::getHardwareCounter(
	const std::string& name, const IHardwareCounter::SpecType& ctrSpec)
{
	// TODO - associate the object we create with the current region
	std::shared_ptr<IHardwareCounter> ret;

	// Check if we have already created a dummy hardware counter set with this
	// name.
	auto iter = allHWCounterSets.find(name);
	if (iter != allHWCounterSets.end()) {
		// We have already created a hw counter set with this name.
		// Return it.
		ret = iter->second;
	}
	else {
		// We have not yet created a hw counter set with this name.
		// Build one and keep track of it.
		// Note with the OSHandler it is always a dummy.
		ret = std::make_shared<dummy::DummyHardwareCounter>(ctrSpec);
		allHWCounterSets[name] = ret;
	}
	return ret;
}

template <typename T, typename V>
void
PerfHandler::collectAllObjectNames(int myRank,
	const std::map<std::string, std::shared_ptr<T>>& myObjs,
	std::map<std::string, PerfObjStatistics<V>>& stats) const
{
	// Collect my own object's names.
	std::vector<std::string> myNames;
	collectMyObjectNames(myObjs, myNames);
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
	auto myNamesBuf = std::make_unique<char[]>(nBytes);
	char* pName = myNamesBuf.get();
	for (auto&& name : myNames) {
		strncpy(pName, name.c_str(), name.length());
		pName += (name.length() + 1); // skip the NUL terminator
	}
	assert(pName == (myNamesBuf.get() + nBytes));

	// Provide all names to root.
	// First, provide the amount of data from each process.
	int cwSize;
	MPI_Comm_size(xolotlComm, &cwSize);
	std::vector<char> allNames((myRank == 0) ? totalNumBytes : 0);
	std::vector<int> allNameCounts((myRank == 0) ? cwSize : 0);
	std::vector<int> allNameDispls((myRank == 0) ? cwSize : 0);

	MPI_Gather(
		&nBytes, 1, MPI_INT, allNameCounts.data(), 1, MPI_INT, 0, xolotlComm);

	// Next, root computes the displacements for data from each process.
	if (myRank == 0) {
		allNameDispls[0] = 0;
		for (int i = 1; i < cwSize; ++i) {
			allNameDispls[i] = allNameDispls[i - 1] + allNameCounts[i - 1];
		}
	}

	// Finally, gather all names to the root process.
	MPI_Gatherv(myNamesBuf.get(), nBytes, MPI_CHAR, allNames.data(),
		allNameCounts.data(), allNameDispls.data(), MPI_CHAR, 0, xolotlComm);

	if (myRank == 0) {
		// Process the gathered names to determine the
		// set of all known object names.
		pName = allNames.data();
		while (pName < (allNames.data() + totalNumBytes)) {
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
		assert(pName == allNames.data() + totalNumBytes);
	}

	// clean up
}

template <typename T>
void
PerfHandler::collectMyObjectNames(
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
PerfHandler::collectMyObjectNames<IHardwareCounter>(
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
PerfHandler::getObjValue(
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
PerfHandler::getObjValue(
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
PerfHandler::aggregateStatistics(int myRank,
	const std::map<std::string, std::shared_ptr<T>>& myObjs,
	std::map<std::string, PerfObjStatistics<V>>& stats) const
{
	// Determine the set of object names known across all processes.
	// Since some processes may define an object that others don't, we
	// have to form the union across all processes.
	// Unfortunately, because the strings are of different lengths,
	// we have a more difficult marshal/unmarshal problem than we'd like.
	collectAllObjectNames<T, V>(myRank, myObjs, stats);
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
		int nameLen = (myRank == 0) ? (int)tsiter->second.name.length() : -1;
		MPI_Bcast(&nameLen, 1, MPI_INT, 0, xolotlComm);
		// we can safely cast away const on the tsiter data string because
		// the only process that accesses that string is rank 0,
		// and it only reads the data.
		std::unique_ptr<char[]> objNamePtr;
		char* objName = nullptr;
		if (myRank == 0) {
			objName = const_cast<char*>(tsiter->second.name.c_str());
		}
		else {
			objNamePtr = std::make_unique<char[]>(nameLen + 1);
			objName = objNamePtr.get();
		}
		MPI_Bcast(objName, nameLen + 1, MPI_CHAR, 0, xolotlComm);

		// do we know about the current object?
		bool knowObject;
		V myVal;
		std::tie<bool, V>(knowObject, myVal) =
			getObjValue<T, V>(myObjs, objName);

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

		// advance to next object
		if (myRank == 0) {
			++tsiter;
		}
	}
}

void
PerfHandler::reportData(std::ostream& os) const
{
	// Output performance information in YAML format.
	os << "\n---\n";

	os << "Rank: " << util::getMPIRank() << '\n';

	const char* indent = "    ";
	os << "\nTimers:\n";
	for (auto&& timer : allTimers) {
		os << indent << timer.first << ": " << timer.second->getValue() << '\n';
	}

	os << "\nCounters:\n";
	for (auto&& ctr : allEventCounters) {
		os << indent << ctr.first << ": " << ctr.second->getValue() << '\n';
	}

	if (allHWCounterSets.empty()) {
		return;
	}
	std::vector<std::string> hwcNames;
	collectMyObjectNames(allHWCounterSets, hwcNames);
	os << "\nHardwareCounters:\n";
	for (auto&& hwc : allHWCounterSets) {
		std::string baseName = hwc.first;

		const auto& spec = hwc.second->getSpecification();
		const auto& vals = hwc.second->getValues();
		for (std::size_t i = 0; i < spec.size(); ++i) {
			auto name = baseName + ":" + hwc.second->getCounterName(spec[i]);
			os << indent << name << ": " << vals[i] << '\n';
		}
	}
}

void
PerfHandler::collectStatistics(PerfObjStatsMap<ITimer::ValType>& timerStats,
	PerfObjStatsMap<IEventCounter::ValType>& counterStats,
	PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats)
{
	// Get the MPI communicator
	auto xolotlComm = util::getMPIComm();
	int myRank;
	MPI_Comm_rank(xolotlComm, &myRank);

	// Aggregate statistics about counters in all processes.
	// First, timers...
	aggregateStatistics<ITimer, ITimer::ValType>(myRank, allTimers, timerStats);

	// ...next event counters...
	aggregateStatistics<IEventCounter, IEventCounter::ValType>(
		myRank, allEventCounters, counterStats);

	// ...finally hardware counters.
	aggregateStatistics<IHardwareCounter, IHardwareCounter::CounterType>(
		myRank, allHWCounterSets, hwCounterStats);
}

void
PerfHandler::reportStatistics(std::ostream& os,
	const PerfObjStatsMap<ITimer::ValType>& timerStats,
	const PerfObjStatsMap<IEventCounter::ValType>& counterStats,
	const PerfObjStatsMap<IHardwareCounter::CounterType>& hwCounterStats) const
{
	// Output performance information in YAML format.
	os << "\n---\n"
	   << "Timers:\n";
	for (auto iter = timerStats.begin(); iter != timerStats.end(); ++iter) {
		iter->second.outputTo(os);
	}

	os << "\nCounters:\n";
	for (auto iter = counterStats.begin(); iter != counterStats.end(); ++iter) {
		iter->second.outputTo(os);
	}

	if (hwCounterStats.empty()) {
		return;
	}
	os << "\nHardwareCounters:\n";
	for (auto iter = hwCounterStats.begin(); iter != hwCounterStats.end();
		 ++iter) {
		iter->second.outputTo(os);
	}
}

void
loadPerfHandlers()
{
}
} // namespace perf
} // namespace xolotl
