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
	 * The bulk temperature in Kelvin
	 */
	double bulkTemperature;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The default constructor is private because the TemperatureHandler
	 * must be initialized with a temperature
	 */
	TemperatureGradientHandler() :
			surfaceTemperature(0.0), bulkTemperature(0.0), dof(0) {
	}

public:

	/**
	 * The constructor
	 *
	 * @param temp The surface temperature
	 * @param grad The bulk temperature
	 */
	TemperatureGradientHandler(double surfTemp, double bulkTemp) :
			surfaceTemperature(surfTemp), bulkTemperature(bulkTemp), dof(0) {
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
	virtual double getTemperature(const NDPoint<3>& fraction, double) const {
		return surfaceTemperature
				+ (bulkTemperature - surfaceTemperature) * fraction[0];
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
	 * This operation sets the heat conductivity to use in the equation.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void setHeatConductivity(double cond) {
		return;
	}

	/**
	 * This operation sets the surface position.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void updateSurfacePosition(int surfacePos) {
		return;
	}

	/**
	 * Compute the flux due to the heat equation.
	 * This method is called by the RHSFunction from the PetscSolver.
	 * Don't do anything.
	 *
	 * \see ITemperatureHandler.h
	 */
	virtual void computeTemperature(double **concVector,
			double *updatedConcOffset, double hxLeft, double hxRight, int xi,
			double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) {
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
			double hxLeft, double hxRight, int xi, double sy = 0.0, int iy = 0,
			double sz = 0.0, int iz = 0) {
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
