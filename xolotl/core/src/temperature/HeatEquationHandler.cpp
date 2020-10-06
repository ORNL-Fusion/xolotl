#include <xolotl/core/Constants.h>
#include <xolotl/core/temperature/HeatEquationHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
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
	double flux, double bulkTemp, int dim) :
	heatFlux(flux),
	bulkTemperature(bulkTemp),
	localTemperature(0.0),
	surfacePosition(0.0),
	heatCoef(0.0),
	heatConductivity(0.0),
	zeroFlux(util::equal(heatFlux, 0.0)),
	dimension(dim),
	oldConcBox(dimension, {0.0, 0.0})
{
}

HeatEquationHandler::HeatEquationHandler(const options::IOptions& options) :
	HeatEquationHandler(options.getTempParam(0), options.getTempParam(1),
		options.getDimensionNumber())
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
		// TODO: Should we throw?
	}
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
HeatEquationHandler::computeTemperature(double** concVector,
	double* updatedConcOffset, double hxLeft, double hxRight, int xi, double sy,
	int iy, double sz, int iz)
{
	if (zeroFlux) {
		return;
	}

	// Initial declaration
	int index = this->_dof;

	// Get the initial concentrations
	double oldConc = concVector[0][index];
	for (int d = 0; d < dimension; ++d) {
		oldConcBox[d][0] = concVector[2 * d + 1][index];
		oldConcBox[d][1] = concVector[2 * d + 2][index];
	}

	double s[3] = {0, sy, sz};

	if (xi == surfacePosition) {
		// Boundary condition with heat flux
		updatedConcOffset[index] += heatCoef * (2.0 / hxLeft) *
			((heatFlux / heatConductivity) +
				(oldConcBox[0][1] - oldConc) / hxRight);
	}
	else {
		// Use a simple midpoint stencil to compute the concentration
		updatedConcOffset[index] += heatCoef * (2.0 / hxLeft) *
			(oldConcBox[0][0] + (hxLeft / hxRight) * oldConcBox[0][1] -
				(1.0 + (hxLeft / hxRight)) * oldConc) /
			(hxLeft + hxRight);
	}

	for (int d = 1; d < dimension; ++d) {
		updatedConcOffset[index] += heatCoef * s[d] *
			(oldConcBox[d][0] + oldConcBox[d][1] - 2.0 * oldConc);
	}
}

bool
HeatEquationHandler::computePartialsForTemperature(double* val, int* indices,
	double hxLeft, double hxRight, int xi, double sy, int iy, double sz, int iz)
{
	if (zeroFlux) {
		return false;
	}

	indices[0] = this->_dof;

	double s[3] = {0, sy, sz};

	val[0] = 1.0 / (hxLeft * hxRight);
	val[1] = 2.0 * heatCoef / (hxLeft * (hxLeft + hxRight));
	val[2] = 2.0 * heatCoef / (hxRight * (hxLeft + hxRight));

	for (int d = 1; d < dimension; ++d) {
		val[0] += s[d];
		val[2 * d + 1] = heatCoef * s[d];
		val[2 * d + 2] = heatCoef * s[d];
	}

	val[0] *= -2.0 * heatCoef;

	if (xi == surfacePosition) {
		val[1] = 0.0;
		val[2] = 2.0 * heatCoef / (hxLeft * hxRight);
	}

	return true;
}
} // namespace temperature
} // namespace core
} // namespace xolotl
