#ifndef HEATEQUATIONHANDLER_H
#define HEATEQUATIONHANDLER_H

#include "ITemperatureHandler.h"
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of the heat equation.
 */
class HeatEquationHandler: public ITemperatureHandler {

private:

	/**
	 * The surface temperature in Kelvin
	 */
	double surfaceTemperature;

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
	double surfacePosition;

	/**
	 * The heat coefficient
	 */
	double heatCoef;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	HeatEquationHandler() :
			surfaceTemperature(0.0), bulkTemperature(0.0), localTemperature(
					0.0), dof(0), surfacePosition(0.0), heatCoef(0.0) {
	}

public:

	/**
	 * The constructor
	 *
	 * @param surfTemp the temperature at the surface
	 * @param bulkTemp the temperature in the bulk
	 */
	HeatEquationHandler(double surfTemp, double bulkTemp) :
			surfaceTemperature(surfTemp), bulkTemperature(bulkTemp), localTemperature(
					0.0), dof(0), surfacePosition(0.0), heatCoef(0.0) {
	}

	/**
	 * The destructor.
	 */
	virtual ~HeatEquationHandler() {
	}

	/**
	 * This operation initializes the ofill and dfill arrays so that the
	 * temperature is connected correctly in the solver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void initializeTemperature(IReactionNetwork *network, int *ofill,
			int *dfill) {
		// Set dof
		dof = network->getDOF();

		// Add the temperature to ofill
		ofill[(dof - 1) * dof + (dof - 1)] = 1;

		// Add the temperature to dfill
		dfill[(dof - 1) * dof + (dof - 1)] = 1;

		return;
	}

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual double getTemperature(const std::vector<double>& position,
			double time) const {
		return xolotlCore::equal(time, 0.0)
				* ((position[0] - surfacePosition < 0.001) * surfaceTemperature
						+ (position[0] - surfacePosition > 0.001)
								* bulkTemperature)
				+ !xolotlCore::equal(time, 0.0) * localTemperature;
	}

	/**
	 * This operation sets the temperature given by the solver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setTemperature(double * solution) {
		localTemperature = solution[dof - 1];
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
	 * This operation sets the surface position.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void updateSurfacePosition(double surfacePos) {
		surfacePosition = surfacePos;
	}

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computeTemperature(double **concVector,
			double *updatedConcOffset, double hxLeft, double hxRight) {
		// Initial declaration
		int index = dof - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index];
		double oldLeftConc = concVector[1][index];
		double oldRightConc = concVector[2][index];

		// Use a simple midpoint stencil to compute the concentration
		double conc = heatCoef * 2.0
				* (oldLeftConc + (hxLeft / hxRight) * oldRightConc
						- (1.0 + (hxLeft / hxRight)) * oldConc)
				/ (hxLeft * (hxLeft + hxRight));

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;

		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computePartialsForTemperature(double *val, int *indices,
			double hxLeft, double hxRight) {
		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[0] = dof - 1;

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, and right grid point
		val[0] = -2.0 * heatCoef / (hxLeft * hxRight); // middle
		val[1] = heatCoef * 2.0 / (hxLeft * (hxLeft + hxRight)); // left
		val[2] = heatCoef * 2.0 / (hxRight * (hxLeft + hxRight)); // right

		return;
	}

};
//end class HeatEquationHandler

}

#endif
