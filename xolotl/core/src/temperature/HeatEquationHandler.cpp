#include <xolotl/core/Constants.h>
#include <xolotl/core/temperature/HeatEquationHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
namespace detail
{
auto heatEqTemperatureHandlerRegistration =
	xolotl::factory::temperature::TemperatureHandlerFactory::Registration<
		HeatEquationHandler>("heat");
}

HeatEquationHandler::HeatEquationHandler(
	double heatFlux, double bulkTemp, int dim, std::string filename) :
	bulkTemperature(bulkTemp),
	localTemperature(0.0),
	surfacePosition(0),
	heatCoef(0.0),
	heatConductivity(0.0),
	zeroFlux(util::equal(heatFlux, 0.0)),
	dimension(dim),
	oldConcBox(dimension, {0.0, 0.0}),
	interfaceLoc(0.0),
	fluxFile(filename)
{
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);
	if (procId == 0) {
		if (fluxFile.size() == 0)
			XOLOTL_LOG << "TemperatureHandler: Using the heat equation with "
						  "a flux of: "
					   << heatFlux << " W nm-2, and a bulk temperature of: "
					   << bulkTemperature << " K";
		else
			XOLOTL_LOG << "TemperatureHandler: Using the heat equation with "
						  "a flux file: "
					   << fluxFile
					   << " , and a bulk temperature of: " << bulkTemperature
					   << " K";
	}

	if (fluxFile.size() > 0) {
		// Open file dataFile.dat containing the time and temperature
		std::ifstream inputFile(fluxFile.c_str());
		std::string line;

		// Read the file and store the values in the two vectors
		while (getline(inputFile, line)) {
			if (!line.length() || line[0] == '#')
				continue;
			double xtemp = 0.0, ytemp = 0.0;
			sscanf(line.c_str(), "%lf %lf", &xtemp, &ytemp);
			time.push_back(xtemp);
			flux.push_back(ytemp);
		}
	}
	else {
		time.push_back(-1.0);
		flux.push_back(heatFlux);
	}
}

HeatEquationHandler::HeatEquationHandler(const options::IOptions& options) :
	HeatEquationHandler(options.getTempParam(0), options.getTempParam(1),
		options.getDimensionNumber(), options.getTempProfileFilename())
{
	// Set the heat coefficient which depends on the material
	auto problemType = options.getMaterial();
	// PSI case
	if (problemType == "W100" || problemType == "W110" ||
		problemType == "W111" || problemType == "W211" ||
		problemType == "TRIDYN") {
		this->setHeatCoefficient(core::tungstenHeatCoefficient);
		this->setHeatConductivity(core::tungstenHeatConductivity);
	}
	// NE case
	else if (problemType == "Fuel") {
		this->setHeatCoefficient(core::uo2HeatCoefficient);
		this->setHeatConductivity(core::uo2HeatConductivity);
	}
	// Fe case
	else if (problemType == "Fe") {
		this->setHeatCoefficient(core::feHeatCoefficient);
		this->setHeatConductivity(core::feHeatConductivity);
	}
	else {
		throw std::runtime_error("\nThe requested material: " + problemType +
			" does not have heat parameters defined for it, cannot use the "
			"temperature option!");
	}

	interfaceLoc = options.getInterfaceLocation();
}

HeatEquationHandler::~HeatEquationHandler()
{
}

double
HeatEquationHandler::getTemperature(
	const plsm::SpaceVector<double, 3>&, double time) const
{
	if (zeroFlux) {
		return bulkTemperature;
	}
	return util::equal(time, 0.0) * bulkTemperature +
		!util::equal(time, 0.0) * localTemperature;
}

void
HeatEquationHandler::computeTemperature(double currentTime, double** concVector,
	double* updatedConcOffset, double hxLeft, double hxRight, int xi, double sy,
	int iy, double sz, int iz)
{
	// Skip if the flux is 0
	if (zeroFlux) {
		return;
	}

	// Initial declaration
	int index = this->_dof;

	auto heatFlux = getHeatFlux(currentTime);

	// Get the initial concentrations
	double oldConc = concVector[0][index];
	for (int d = 0; d < dimension; ++d) {
		oldConcBox[d][0] = concVector[2 * d + 1][index];
		oldConcBox[d][1] = concVector[2 * d + 2][index];
	}

	// Adjust the parameters
	double alpha = getLocalHeatAlpha(xi), beta = getLocalHeatBeta(oldConc),
		   gamma = getLocalHeatGamma(oldConc), dAlpha = getDAlpha(xi),
		   dBeta = getDBeta(oldConc);

	double s[3] = {0, sy, sz};

	double x1 = xGrid[xi] - xGrid[surfacePosition + 1];
	double x2 = xGrid[xi + 1] - xGrid[surfacePosition + 1];

	// Surface and interface
	if (xi == surfacePosition or (interfaceLoc > x1 and interfaceLoc <= x2)) {
		// Boundary condition with heat flux
		updatedConcOffset[index] += (2.0 * heatFlux * gamma / hxLeft) +
			(2.0 * alpha * beta * gamma) * (oldConcBox[0][1] - oldConc) /
				(hxLeft * hxRight);
		// Second term for temperature dependent conductivity
		updatedConcOffset[index] += -heatFlux * dAlpha * gamma / alpha +
			heatFlux * heatFlux * gamma * dBeta / (alpha * beta * beta);
	}
	else if (xi == bulkPosition) {
		// Boundary condition with heat flux
		double bulkHeatFlux = getBulkHeatFlux(oldConc);
		updatedConcOffset[index] += -(2.0 * bulkHeatFlux * gamma / hxRight) +
			(2.0 * alpha * beta * gamma) * (oldConcBox[0][0] - oldConc) /
				(hxLeft * hxRight);
		// Second term for temperature dependent conductivity
		updatedConcOffset[index] += -bulkHeatFlux * dAlpha * gamma / alpha +
			bulkHeatFlux * bulkHeatFlux * gamma * dBeta / (alpha * beta * beta);
	}
	else {
		// Use a simple midpoint stencil to compute the concentration
		updatedConcOffset[index] += 2.0 * alpha * beta * gamma *
			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
				(1.0 + (hxLeft / hxRight)) * oldConc) /
			(hxLeft * (hxLeft + hxRight));
		// Second term for temperature dependent conductivity
		updatedConcOffset[index] += dAlpha * beta * gamma *
				(oldConcBox[0][1] - oldConcBox[0][0]) / (hxLeft + hxRight) +
			alpha * dBeta * gamma * (oldConcBox[0][1] - oldConcBox[0][0]) *
				(oldConcBox[0][1] - oldConcBox[0][0]) /
				((hxLeft + hxRight) * (hxLeft + hxRight));
	}

	// Deal with the potential additional dimensions
	for (int d = 1; d < dimension; ++d) {
		updatedConcOffset[index] += alpha * beta * gamma * s[d] *
			(oldConcBox[d][0] + oldConcBox[d][1] - 2.0 * oldConc);
	}
}

bool
HeatEquationHandler::computePartialsForTemperature(double currentTime,
	double** concVector, double* val, IdType* indices, double hxLeft,
	double hxRight, int xi, double sy, int iy, double sz, int iz)
{
	// Skip if the flux is 0
	if (zeroFlux) {
		return false;
	}

	// Initial declaration
	int index = this->_dof;

	auto heatFlux = getHeatFlux(currentTime);

	// Get the initial concentrations
	double oldConc = concVector[0][index];
	for (int d = 0; d < dimension; ++d) {
		oldConcBox[d][0] = concVector[2 * d + 1][index];
		oldConcBox[d][1] = concVector[2 * d + 2][index];
	}

	// Get the DOF
	indices[0] = index;

	double s[3] = {0, sy, sz};

	double alpha = getLocalHeatAlpha(xi), beta = getLocalHeatBeta(oldConc),
		   gamma = getLocalHeatGamma(oldConc), dAlpha = getDAlpha(xi),
		   dBeta = getDBeta(oldConc), dGamma = getDGamma(oldConc),
		   ddBeta = getDDBeta(oldConc);

	// Compute the partials along the depth
	val[0] = -2.0 * alpha * beta * gamma / (hxLeft * hxRight) +
		alpha * (ddBeta * gamma + dGamma * dBeta) *
			(oldConcBox[0][1] - oldConcBox[0][0]) *
			(oldConcBox[0][1] - oldConcBox[0][0]) /
			((hxLeft + hxRight) * (hxLeft + hxRight)) +
		2.0 * alpha * (dBeta * gamma + beta * dGamma) *
			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
				(1.0 + (hxLeft / hxRight)) * oldConc) /
			((hxLeft + hxRight) * hxLeft) +
		dAlpha * (dBeta * gamma + dGamma * beta) *
			(oldConcBox[0][1] - oldConcBox[0][0]) / (hxLeft + hxRight);
	val[1] = 2.0 * alpha * beta * gamma / (hxLeft * (hxLeft + hxRight)) -
		beta * dAlpha * gamma / (hxLeft + hxRight) +
		2.0 * alpha * dBeta * gamma * (oldConcBox[0][0] - oldConcBox[0][1]) /
			((hxLeft + hxRight) * (hxLeft + hxRight));
	val[2] = 2.0 * alpha * beta * gamma / (hxRight * (hxLeft + hxRight)) +
		beta * dAlpha * gamma / (hxLeft + hxRight) +
		2.0 * alpha * dBeta * gamma * (oldConcBox[0][1] - oldConcBox[0][0]) /
			((hxLeft + hxRight) * (hxLeft + hxRight));

	// Deal with the potential additional dimensions
	for (int d = 1; d < dimension; ++d) {
		val[0] -= 2.0 * alpha * beta * gamma * s[d];
		val[2 * d + 1] = alpha * beta * gamma * s[d];
		val[2 * d + 2] = alpha * beta * gamma * s[d];
	}

	double x1 = xGrid[xi] - xGrid[surfacePosition + 1];
	double x2 = xGrid[xi + 1] - xGrid[surfacePosition + 1];

	// Boundary condition with the heat flux
	if (xi == surfacePosition or (interfaceLoc > x1 and interfaceLoc <= x2)) {
		val[0] = 2.0 * heatFlux * dGamma / hxLeft -
			2.0 * alpha * beta * gamma / (hxLeft * hxRight) +
			2.0 * alpha * (dBeta * gamma + dGamma * beta) *
				(oldConcBox[0][1] - oldConc) / (hxLeft * hxRight) -
			heatFlux * dGamma * dAlpha / alpha +
			heatFlux * heatFlux * (gamma * ddBeta + dGamma * dBeta) /
				(alpha * beta * beta) -
			2.0 * heatFlux * heatFlux * gamma * dBeta * dBeta *
				(alpha * beta * beta * beta);
		val[1] = 0.0;
		val[2] = 2.0 * alpha * beta * gamma / (hxLeft * hxRight);
	}
	else if (xi == bulkPosition) {
		double bulkHeatFlux = getBulkHeatFlux(oldConc);
		double dBulk = getBulkHeatFluxDerivative(oldConc);
		val[0] = -2.0 * bulkHeatFlux * dGamma / hxRight -
			2.0 * alpha * beta * gamma / (hxLeft * hxRight) +
			2.0 * alpha * (dBeta * gamma + dGamma * beta) *
				(oldConcBox[0][0] - oldConc) / (hxLeft * hxRight) -
			bulkHeatFlux * dAlpha * dGamma / alpha +
			bulkHeatFlux * bulkHeatFlux * (gamma * ddBeta + dGamma * dBeta) /
				(alpha * beta * beta) -
			2.0 * bulkHeatFlux * bulkHeatFlux * gamma * dBeta * dBeta *
				(alpha * beta * beta * beta) -
			dBulk * gamma * dAlpha / alpha - 2.0 * dBulk * gamma / hxRight +
			2.0 * bulkHeatFlux * dBulk * gamma * dBeta / (alpha * beta * beta);
		val[1] = 2.0 * alpha * beta * gamma / (hxLeft * hxRight);
		val[2] = 0.0;
	}

	return true;
}

double
HeatEquationHandler::getLocalHeatAlpha(int xi) const
{
	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	if (x < interfaceLoc)
		return 0.2;
	else
		return 1.0;

	double heatCond = (0.2 + 0.8 / (1.0 + exp((interfaceLoc - x) * 0.1)));
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatBeta(double temp) const
{
	double heatCond =
		(-2.4e-12 * temp * temp * temp * temp - 3.67e-10 * temp * temp * temp +
			4.65e-5 * temp * temp - 0.131 * temp + 207.0) *
		1.0e-9;
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatGamma(double temp) const
{
	double tungstenDensity = 19.25e-24; // kg nm-3
	double capacity = -4.01e-6 * temp * temp + 3.61e-2 * temp + 119.0;
	return 1.0 / (tungstenDensity * capacity);
}

double
HeatEquationHandler::getDBeta(double temp) const
{
	double heatCond = (-9.6e-12 * temp * temp * temp - 11.01e-10 * temp * temp +
						  9.3e-5 * temp - 0.131) *
		1.0e-9;
	return heatCond;
}

double
HeatEquationHandler::getDDBeta(double temp) const
{
	double heatCond =
		(-28.8e-12 * temp * temp - 22.02e-10 * temp + 9.3e-5) * 1.0e-9;
	return heatCond;
}

double
HeatEquationHandler::getDAlpha(int xi) const
{
	return 0.0;
	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	double heatCond = 0.1 * 0.8 * exp((interfaceLoc - x) * 0.1) /
		((1.0 + exp((interfaceLoc - x) * 0.1)) *
			(1.0 + exp((interfaceLoc - x) * 0.1)));
	if (interfaceLoc - x > 1000.0)
		return 0.0;
	return heatCond;
}

double
HeatEquationHandler::getDGamma(double temp) const
{
	double tungstenDensity = 19.25e-24; // kg nm-3
	double capacity = -4.01e-6 * temp * temp + 3.61e-2 * temp + 119.0;
	double dCapacity = 3.61e-2 - 8.02e-6 * temp;
	return -dCapacity / (tungstenDensity * capacity * capacity);
}

double
HeatEquationHandler::getHeatFlux(double currentTime)
{
	// Initialize the value to return
	double Qinter = 6.0e-12;
	double Q0 = 7.4e-10;
	double tau = 0.8 * 2.5e-4;
	double freq = 50.0;
	double cycleTime =
		currentTime - (std::floor(currentTime * freq)) * (1.0 / freq);
	double x = 0.0;
	if (cycleTime > 0.0)
		x = tau / cycleTime;
	//	return Qinter;
	return Qinter + Q0 * (1.0 + x * x) * x * x * exp(-x * x);

	double f = 0.0;
	// If the time is smaller than or equal than the first stored time
	if (currentTime <= time[0])
		return f = flux[0];

	// If the time is larger or equal to the last stored time
	if (currentTime >= time[time.size() - 1])
		return f = flux[time.size() - 1];

	// Else loop to determine the interval the time falls in
	// i.e. time[k] < time < time[k + 1]
	for (unsigned int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k])
			continue;
		if (currentTime > time[k + 1])
			continue;

		// Compute the amplitude following a linear interpolation between
		// the two stored values
		f = flux[k] +
			(flux[k + 1] - flux[k]) * (currentTime - time[k]) /
				(time[k + 1] - time[k]);
		break;
	}

	return f;
}

double
HeatEquationHandler::getBulkHeatFlux(double temp) const
{
	// Convert the temperature to Celsius
	double tempCelsius = temp - 273.0;

	return (-8.7993 + 0.0845 * tempCelsius) * 1.0e-12;
}

double
HeatEquationHandler::getBulkHeatFluxDerivative(double temp) const
{
	return 0.0845 * 1.0e-12;
}

} // namespace temperature
} // namespace core
} // namespace xolotl
