#include <fstream>

#include <xolotl/core/temperature/ProfileHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
namespace detail
{
auto profileTemperatureHandlerRegistration =
	xolotl::factory::temperature::TemperatureHandlerFactory::Registration<
		ProfileHandler>("profile");
}

ProfileHandler::ProfileHandler(const std::string& profileFileName) :
	tempFile(profileFileName)
{
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);
	if (procId == 0) {
		XOLOTL_LOG << "TemperatureHandler: Using the time profile defined in: "
				   << tempFile;
	}
}

ProfileHandler::ProfileHandler(const options::IOptions& options) :
	ProfileHandler(options.getTempProfileFilename())
{
}

ProfileHandler::~ProfileHandler()
{
}

void
ProfileHandler::initialize(const int dof)
{
	TemperatureHandler::initialize(dof);

	// Open file dataFile.dat containing the time and temperature
	std::ifstream inputFile(tempFile.c_str());
	std::string line;

	// Read the file and store the values in the two vectors
	while (getline(inputFile, line)) {
		if (!line.length() || line[0] == '#')
			continue;
		double xtemp = 0.0, ytemp = 0.0;
		sscanf(line.c_str(), "%lf %lf", &xtemp, &ytemp);
		time.push_back(xtemp);
		temp.push_back(ytemp);
	}
}

double
ProfileHandler::getTemperature(
	const plsm::SpaceVector<double, 3>&, double currentTime) const
{
	// Initialize the value to return
	double f = 0.0;

	// If the time is smaller than or equal than the first stored time
	if (currentTime <= time[0])
		return temp[0];

	// If the time is larger or equal to the last stored time
	if (currentTime >= time[time.size() - 1])
		return temp[time.size() - 1];

	// Else loop to determine the interval the time falls in
	// i.e. time[k] < time < time[k + 1]
	for (unsigned int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k])
			continue;
		if (currentTime > time[k + 1])
			continue;

		// Compute the amplitude following a linear interpolation between
		// the two stored values
		f = temp[k] +
			(temp[k + 1] - temp[k]) * (currentTime - time[k]) /
				(time[k + 1] - time[k]);
		break;
	}

	return f;
}
} // namespace temperature
} // namespace core
} // namespace xolotl
