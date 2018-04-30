#ifndef TEMPERATUREGRADIENTHANDLER_H
#define TEMPERATUREGRADIENTHANDLER_H

#include "ITemperatureHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITemperatureHandler, it is responsible for the
 * handling of a temperature constant with time but changing with location.
 */
class TemperatureGradientHandler: public ITemperatureHandler {

private:

	/**
	 * The surface temperature in Kelvin
	 */
	double surfaceTemperature;

	/**
	 * The temperature gradient
	 */
	double gradient;

	/**
	 * The surface position
	 */
	double surfacePosition;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	TemperatureGradientHandler() :
			surfaceTemperature(0.0), gradient(0.0), surfacePosition(0.0), dof(0) {
	}

public:

	/**
	 * The constructor
	 *
	 * @param temp The surface temperature
	 * @param grad The temperature gradient
	 */
	TemperatureGradientHandler(double temp, double grad) :
			surfaceTemperature(temp), gradient(grad), surfacePosition(0.0), dof(
					0) {
	}

	/**
	 * The destructor.
	 */
	virtual ~TemperatureGradientHandler() {
	}

	/**
	 * This operation initializes the ofill and dfill arrays so that the
	 * temperature is connected correctly in the solver.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void initializeTemperature(const IReactionNetwork& network,
            IReactionNetwork::SparseFillMap& ofillMap,
            IReactionNetwork::SparseFillMap& dfillMap) {

		// Set dof
		dof = network.getDOF();

		// Add the temperature to ofill
        ofillMap[(dof - 1)].emplace_back(dof - 1);

		// Add the temperature to dfill
        dfillMap[(dof - 1)].emplace_back(dof - 1);

		return;
	}

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 *
	 * @return The temperature
	 */
	virtual double getTemperature(const Point<3>& position, double) const {
		return surfaceTemperature - (position[0] - surfacePosition) * gradient;
	}

	/**
	 * This operation sets the temperature given by the solver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setTemperature(double * solution) {
		return;
	}

	/**
	 * This operation sets the heat coefficient to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setHeatCoefficient(double coef) {
		return;
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
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computeTemperature(double **concVector,
			double *updatedConcOffset, double hxLeft, double hxRight) {
		return;
	}

	/**
	 * Compute the partials due to the heat equation.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * Don't do anything.
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
		val[0] = 0.0; // middle
		val[1] = 0.0; // left
		val[2] = 0.0; // right

		return;
	}

};
//end class TemperatureGradientHandler

}

#endif
