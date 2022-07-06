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
	double midHeatCoef = getLocalHeatCoefficient(xi, oldConc),
		   midHeatCond = getLocalHeatConductivity(xi, oldConc);
	double alpha = getLocalHeatAlpha(xi), beta = getLocalHeatBeta(oldConc),
		   dAlpha = getLocalHeatCondSpatialDerivative(xi),
		   dBeta = getLocalHeatCondTempDerivative(oldConc);

	double s[3] = {0, sy, sz};

	// Surface and interface
	//	if (xi == surfacePosition ||
	//		fabs(xGrid[xi + 1] - xGrid[surfacePosition + 1] - interfaceLoc)
	//< 2.0)
	//{
	//		// Boundary condition with heat flux
	//		updatedConcOffset[index] += midHeatCoef * (2.0 / hxLeft) *
	//			((heatFlux / midHeatCond) + (oldConcBox[0][1] - oldConc) /
	// hxRight);
	//		// Second term for temperature dependent conductivity
	//		updatedConcOffset[index] += midHeatCoef * heatFlux * heatFlux *
	//			getLocalHeatCondTempDerivative(xi, oldConc) /
	//			(midHeatCond * midHeatCond * midHeatCond);
	//	}
	//	else {
	//		// Use a simple midpoint stencil to compute the concentration
	//		updatedConcOffset[index] += midHeatCoef * (2.0 / hxLeft) *
	//			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
	//				(1.0 + (hxLeft / hxRight)) * oldConc) /
	//			(hxLeft + hxRight);
	//		// Second term for temperature dependent conductivity
	//		updatedConcOffset[index] += midHeatCoef *
	//			(oldConcBox[0][1] - oldConcBox[0][0]) *
	//			(oldConcBox[0][1] - oldConcBox[0][0]) *
	//			getLocalHeatCondTempDerivative(xi, oldConc) /
	//			(midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight));
	//	}
	if (xi == surfacePosition) {
		// Boundary condition with heat flux
		updatedConcOffset[index] += midHeatCoef * (2.0 / hxLeft) *
			((heatFlux / midHeatCond) + (oldConcBox[0][1] - oldConc) / hxRight);
		// Second term for temperature dependent conductivity
		updatedConcOffset[index] += -midHeatCoef * heatFlux * beta * dAlpha /
				(midHeatCond * midHeatCond) +
			midHeatCoef * heatFlux * heatFlux * alpha * dBeta /
				(midHeatCond * midHeatCond * midHeatCond);
	}
	else {
		// Use a simple midpoint stencil to compute the concentration
		updatedConcOffset[index] += midHeatCoef * (2.0 / hxLeft) *
			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
				(1.0 + (hxLeft / hxRight)) * oldConc) /
			(hxLeft + hxRight);
		// Second term for temperature dependent conductivity
		updatedConcOffset[index] += midHeatCoef *
				(oldConcBox[0][1] - oldConcBox[0][0]) * beta * dAlpha /
				(midHeatCond * (hxLeft + hxRight)) +
			midHeatCoef * alpha * dBeta *
				(oldConcBox[0][1] - oldConcBox[0][0]) *
				(oldConcBox[0][1] - oldConcBox[0][0]) /
				(midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight));
	}

	// Deal with the potential additional dimensions
	for (int d = 1; d < dimension; ++d) {
		updatedConcOffset[index] += midHeatCoef * s[d] *
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

	double midHeatCoef = getLocalHeatCoefficient(xi, oldConc),
		   midHeatCond = getLocalHeatConductivity(xi, oldConc);
	double alpha = getLocalHeatAlpha(xi), beta = getLocalHeatBeta(oldConc),
		   dAlpha = getLocalHeatCondSpatialDerivative(xi),
		   dBeta = getLocalHeatCondTempDerivative(oldConc),
		   ddBeta = getLocalHeatCondTempSecondDerivative(oldConc);

	// Compute the partials along the depth
	//	val[0] = -2.0 * midHeatCoef / (hxLeft * hxRight) +
	//		midHeatCoef * getLocalHeatCondTempSecondDerivative(xi, oldConc) *
	//			(oldConcBox[0][1] - oldConcBox[0][0]) *
	//			(oldConcBox[0][1] - oldConcBox[0][0]) /
	//			(midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight)) +
	//		midHeatCoef * (2.0 / hxLeft) *
	//			getLocalHeatCondTempDerivative(xi, oldConc) *
	//			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
	//				(1.0 + (hxLeft / hxRight)) * oldConc) /
	//			((hxLeft + hxRight) * midHeatCond);
	//	val[1] = 2.0 * midHeatCoef / (hxLeft * (hxLeft + hxRight)) +
	//		2.0 * midHeatCoef * (oldConcBox[0][0] - oldConcBox[0][1]) *
	//			getLocalHeatCondTempDerivative(xi, oldConc) /
	//			(midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight));
	//	val[2] = 2.0 * midHeatCoef / (hxRight * (hxLeft + hxRight)) +
	//		2.0 * midHeatCoef * (oldConcBox[0][1] - oldConcBox[0][0]) *
	//			getLocalHeatCondTempDerivative(xi, oldConc) /
	//			(midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight));
	val[0] = -2.0 * midHeatCoef / (hxLeft * hxRight) +
		midHeatCoef * alpha * ddBeta * (oldConcBox[0][1] - oldConcBox[0][0]) *
			(oldConcBox[0][1] - oldConcBox[0][0]) /
			(midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight)) +
		midHeatCoef * (2.0 / hxLeft) * alpha * dBeta *
			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
				(1.0 + (hxLeft / hxRight)) * oldConc) /
			((hxLeft + hxRight) * midHeatCond) +
		midHeatCoef * dBeta * dAlpha * (oldConcBox[0][1] - oldConcBox[0][0]) /
			((hxLeft + hxRight) * midHeatCond);
	val[1] = 2.0 * midHeatCoef / (hxLeft * (hxLeft + hxRight)) -
		midHeatCoef * beta * dAlpha / (midHeatCond * (hxLeft + hxRight)) +
		2.0 * midHeatCoef * (oldConcBox[0][0] - oldConcBox[0][1]) * alpha *
			dBeta / (midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight));
	val[2] = 2.0 * midHeatCoef / (hxRight * (hxLeft + hxRight)) +
		midHeatCoef * beta * dAlpha / (midHeatCond * (hxLeft + hxRight)) +
		2.0 * midHeatCoef * (oldConcBox[0][1] - oldConcBox[0][0]) * alpha *
			dBeta / (midHeatCond * (hxLeft + hxRight) * (hxLeft + hxRight));

	// Deal with the potential additional dimensions
	for (int d = 1; d < dimension; ++d) {
		val[0] -= 2.0 * midHeatCoef * s[d];
		val[2 * d + 1] = midHeatCoef * s[d];
		val[2 * d + 2] = midHeatCoef * s[d];
	}

	// Boundary condition with the heat flux
	//	if (xi == surfacePosition ||
	//		fabs(xGrid[xi + 1] - xGrid[surfacePosition + 1] - interfaceLoc)
	//< 2.0) { 		val[0] = -2.0 * midHeatCoef / (hxLeft * hxRight)
	//+ 			2.0 * midHeatCoef * getLocalHeatCondTempDerivative(xi,
	// oldConc) * 				(oldConcBox[0][1] - oldConc)
	/// 				(midHeatCond * hxLeft * hxRight) + 			heatFlux
	/// * heatFlux * midHeatCoef *
	//				getLocalHeatCondTempSecondDerivative(xi, oldConc) /
	//				(midHeatCond * midHeatCond * midHeatCond) -
	//			2.0 * heatFlux * heatFlux * midHeatCoef *
	//				getLocalHeatCondTempDerivative(xi, oldConc) *
	//				getLocalHeatCondTempDerivative(xi, oldConc) *
	//				(midHeatCond * midHeatCond * midHeatCond * midHeatCond);
	//		val[1] = 0.0;
	//		val[2] = 2.0 * midHeatCoef / (hxLeft * hxRight);
	//	}
	if (xi == surfacePosition) {
		val[0] = -2.0 * midHeatCoef / (hxLeft * hxRight) +
			2.0 * midHeatCoef * alpha * dBeta * (oldConcBox[0][1] - oldConc) /
				(midHeatCond * hxLeft * hxRight) +
			heatFlux * heatFlux * midHeatCoef * alpha * ddBeta /
				(midHeatCond * midHeatCond * midHeatCond) -
			2.0 * heatFlux * heatFlux * midHeatCoef * alpha * alpha * dBeta *
				dBeta *
				(midHeatCond * midHeatCond * midHeatCond * midHeatCond) +
			midHeatCoef * heatFlux * beta * alpha * dAlpha * dBeta /
				(midHeatCond * midHeatCond * midHeatCond) -
			midHeatCoef * heatFlux * dBeta * dAlpha /
				(midHeatCond * midHeatCond);
		val[1] = 0.0;
		val[2] = 2.0 * midHeatCoef / (hxLeft * hxRight);
	}

	return true;
}

double
HeatEquationHandler::getLocalHeatConductivity(int xi, double temp) const
{
	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	double lnT = log(temp);
	double heatCond = (A * lnT * lnT + B * lnT + C) * 1.0e-9 *
		(0.2 + 0.8 / (1.0 + exp((interfaceLoc - x) * 0.1)));
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatAlpha(int xi) const
{
	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	double heatCond = (0.2 + 0.8 / (1.0 + exp((interfaceLoc - x) * 0.1)));
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatBeta(double temp) const
{
	double lnT = log(temp);
	double heatCond = (A * lnT * lnT + B * lnT + C) * 1.0e-9;
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatCondTempDerivative(double temp) const
{
	double lnT = log(temp);
	double heatCond = (B + 2.0 * A * lnT) * 1.0e-9 / temp;
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatCondTempSecondDerivative(double temp) const
{
	double lnT = log(temp);
	double heatCond = (-B + 2.0 * A * (1.0 - lnT)) * 1.0e-9 / (temp * temp);
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatCondSpatialDerivative(int xi) const
{
	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	double heatCond = 0.1 * 0.8 * exp((interfaceLoc - x) * 0.1) /
		((1.0 + exp((interfaceLoc - x) * 0.1)) *
			(1.0 + exp((interfaceLoc - x) * 0.1)));
	if (interfaceLoc - x > 1000.0)
		return 0.0;
	return heatCond;
}

double
HeatEquationHandler::getLocalHeatCoefficient(int xi, double temp) const
{
	return heatCoef * getLocalHeatConductivity(xi, temp);
}

double
HeatEquationHandler::getHeatFlux(double currentTime)
{
	// Initialize the value to return
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
} // namespace temperature
} // namespace core
} // namespace xolotl
