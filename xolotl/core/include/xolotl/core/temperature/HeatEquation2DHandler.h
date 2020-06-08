#ifndef HEATEQUATION2DHANDLER_H
#define HEATEQUATION2DHANDLER_H

#include <xolotl/core/temperature/ITemperatureHandler.h>
#include <xolotl/util/MathUtils.h>
#include <xolotl/core/Constants.h>

namespace xolotl {
namespace core {
namespace temperature {

/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of the heat equation in 2D.
 */
class HeatEquation2DHandler: public ITemperatureHandler {

private:

	/**
	 * The heat flux in W.m-2
	 */
	double heatFlux;

	/**
	 * The bulk temperature in Kelvin
	 */
	double bulkTemperature;

	/**
	 * The local temperature in Kelvin
	 */
	double localTemperature;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The surface position
	 */
	int surfacePosition;

	/**
	 * The heat coefficient
	 */
	double heatCoef;

	/**
	 * The heat conductivity
	 */
	double heatConductivity;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	HeatEquation2DHandler() :
			heatFlux(0.0), bulkTemperature(0.0), localTemperature(0.0), dof(0), surfacePosition(
					0.0), heatCoef(0.0), heatConductivity(0.0) {
	}

public:

	/**
	 * The constructor
	 *
	 * @param flux The heat flux
	 * @param bulkTemp The temperature in the bulk
	 */
	HeatEquation2DHandler(double flux, double bulkTemp) :
			heatFlux(flux), bulkTemperature(bulkTemp), localTemperature(0.0), dof(
					0), surfacePosition(0.0), heatCoef(0.0), heatConductivity(
					0.0) {
	}

	/**
	 * The destructor.
	 */
	virtual ~HeatEquation2DHandler() {
	}

	/**
	 * This operation initializes the ofill and dfill arrays so that the
	 * temperature is connected correctly in the solver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void initializeTemperature(const int _dof,
			network::IReactionNetwork::SparseFillMap& ofillMap,
			network::IReactionNetwork::SparseFillMap& dfillMap) {

		// Set dof
		dof = _dof;

		// Add the temperature to ofill
		ofillMap[dof].emplace_back(dof);

		// Add the temperature to dfill
		dfillMap[dof].emplace_back(dof);

		return;
	}

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual double getTemperature(const plsm::SpaceVector<double, 3>&, double time) const {
		return util::equal(time, 0.0) * bulkTemperature
				+ !util::equal(time, 0.0) * localTemperature;
	}

	/**
	 * This operation sets the temperature given by the solver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setTemperature(double * solution) {
		localTemperature = solution[dof];
	}

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setHeatCoefficient(double coef) {
		heatCoef = coef;
	}

	/**
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setHeatConductivity(double cond) {
		heatConductivity = cond;
	}

	/**
	 * This operation sets the surface position.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void updateSurfacePosition(int surfacePos) {
		surfacePosition = surfacePos;
	}

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computeTemperature(double **concVector,
			double *updatedConcOffset, double hxLeft, double hxRight, int xi,
			double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) {
		// Initial declaration
		int index = dof;

		// Get the initial concentrations
		double oldConc = concVector[0][index];
		double oldLeftConc = concVector[1][index];
		double oldRightConc = concVector[2][index];
		double oldBottomConc = concVector[3][index];
		double oldTopConc = concVector[4][index];

		// Boundary condition with heat flux
		if (xi == surfacePosition) {
			// Include the flux boundary condition
			updatedConcOffset[index] +=
					heatCoef
							* ((2.0 / hxLeft)
									* ((heatFlux / heatConductivity)
											+ (oldRightConc - oldConc) / hxRight)
									+ sy
											* (oldBottomConc + oldTopConc
													- 2.0 * oldConc));

			return;
		}

		// Use a simple midpoint stencil to compute the concentration
		updatedConcOffset[index] += heatCoef
				* (2.0
						* (oldLeftConc + (hxLeft / hxRight) * oldRightConc
								- (1.0 + (hxLeft / hxRight)) * oldConc)
						/ (hxLeft * (hxLeft + hxRight))
						+ sy * (oldBottomConc + oldTopConc - 2.0 * oldConc));

		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual bool computePartialsForTemperature(double *val, int *indices,
			double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
			double sz = 0.0, int iz = 0) {
		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[0] = dof;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, and right grid point
		val[0] = -2.0 * heatCoef * (1.0 / (hxLeft * hxRight) + sy); // middle
		val[1] = heatCoef * 2.0 / (hxLeft * (hxLeft + hxRight)); // left
		val[2] = heatCoef * 2.0 / (hxRight * (hxLeft + hxRight)); // right
		val[3] = heatCoef * sy; // bottom
		val[4] = heatCoef * sy; // top

		if (xi == surfacePosition) {
			val[1] = 0.0;
			val[2] = 2.0 * heatCoef / (hxLeft * hxRight);
		}

		return true;
	}

};
//end class HeatEquation2DHandler

}
}
}

#endif
