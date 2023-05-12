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
	double flux, double bulkTemp, int dim) :
	heatFlux(flux),
	bulkTemperature(bulkTemp),
	localTemperature(0.0),
	surfacePosition(0.0),
	heatCoef(0.0),
	heatConductivity(0.0),
	zeroFlux(util::equal(heatFlux, 0.0)),
	dimension(dim)
{
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);
	if (procId == 0) {
		XOLOTL_LOG << "TemperatureHandler: Using the heat equation with "
					  "a flux of: "
				   << heatFlux
				   << " W nm-2, and a bulk temperature of: " << bulkTemperature
				   << " K";
	}
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
		throw std::runtime_error("\nThe requested material: " + problemType +
			" does not have heat parameters defined for it, cannot use the "
			"temperature option!");
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
HeatEquationHandler::setTemperature(Kokkos::View<const double*> solution)
{
	auto elemView = subview(solution, this->_dof);
	auto elem_h = create_mirror_view(elemView);
	deep_copy(elem_h, elemView);

	localTemperature = elem_h();
}

void
HeatEquationHandler::computeTemperature(Kokkos::View<const double*>* concVector,
	Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
	int xi, double sy, int iy, double sz, int iz)
{
	// Skip if the flux is 0
	if (zeroFlux) {
		return;
	}

	// Initial declaration
	int index = this->_dof;
	Kokkos::Array<double, 3> s = {0, sy, sz};
	Kokkos::Array<Kokkos::View<const double*>, 7> concVec;
	concVec[0] = concVector[0];
	for (int d = 0; d < dimension; ++d) {
		concVec[2 * d + 1] = concVector[2 * d + 1];
		concVec[2 * d + 2] = concVector[2 * d + 2];
	}

	auto dim = dimension;
	auto surfPos = surfacePosition;
	auto htCoef = heatCoef;
	auto htFlux = heatFlux;
	auto htConductivity = heatConductivity;

	////////////////////////////////////////////////////////////////////////////
	// TODO: Move this to the constructor (or change member to be a View
	// util::Array<double, 3, 2> oldBox;
	////////////////////////////////////////////////////////////////////////////

	Kokkos::parallel_for(
		1, KOKKOS_LAMBDA(std::size_t) {
			// Get the initial concentrations
			double oldConc = concVec[0][index];
			util::Array<double, 3, 2> oldBox;
			for (int d = 0; d < dim; ++d) {
				oldBox[d][0] = concVec[2 * d + 1][index];
				oldBox[d][1] = concVec[2 * d + 2][index];
			}

			if (xi == surfPos) {
				// Boundary condition with heat flux
				updatedConcOffset[index] += htCoef * (2.0 / hxLeft) *
					((htFlux / htConductivity) +
						(oldBox[0][1] - oldConc) / hxRight);
			}
			else {
				// Use a simple midpoint stencil to compute the concentration
				updatedConcOffset[index] += htCoef * (2.0 / hxLeft) *
					(oldBox[0][0] + (hxLeft / hxRight) * oldBox[0][1] -
						(1.0 + (hxLeft / hxRight)) * oldConc) /
					(hxLeft + hxRight);
			}

			// Deal with the potential additional dimensions
			for (int d = 1; d < dim; ++d) {
				updatedConcOffset[index] += htCoef * s[d] *
					(oldBox[d][0] + oldBox[d][1] - 2.0 * oldConc);
			}
		});
}

bool
HeatEquationHandler::computePartialsForTemperature(double* val, IdType* indices,
	double hxLeft, double hxRight, int xi, double sy, int iy, double sz, int iz)
{
	// Skip if the flux is 0
	if (zeroFlux) {
		return false;
	}

	// Get the DOF
	indices[0] = this->_dof;

	double s[3] = {0, sy, sz};

	// Compute the partials along the depth
	val[0] = 1.0 / (hxLeft * hxRight);
	val[1] = 2.0 * heatCoef / (hxLeft * (hxLeft + hxRight));
	val[2] = 2.0 * heatCoef / (hxRight * (hxLeft + hxRight));

	// Deal with the potential additional dimensions
	for (int d = 1; d < dimension; ++d) {
		val[0] += s[d];
		val[2 * d + 1] = heatCoef * s[d];
		val[2 * d + 2] = heatCoef * s[d];
	}

	val[0] *= -2.0 * heatCoef;

	// Boundary condition with the heat flux
	if (xi == surfacePosition) {
		val[1] = 0.0;
		val[2] = 2.0 * heatCoef / (hxLeft * hxRight);
	}

	return true;
}
} // namespace temperature
} // namespace core
} // namespace xolotl
