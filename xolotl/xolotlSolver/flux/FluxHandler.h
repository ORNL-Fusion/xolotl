#ifndef FLUXHANDLER_H
#define FLUXHANDLER_H

#include "IFluxHandler.h"
#include <vector>

namespace xolotlSolver {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing flux calculations.
 */
class FluxHandler: public IFluxHandler {

private:

protected:

	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position)
	 */
	std::vector<double> incidentFluxVec;

	/**
	 * Step size between each grid point
	 */
	double stepSize;

	/**
	 * Helium fluence
	 */
	double heFluence;

	/**
	 * Should the program use the maximum Helium fluence value?
	 */
	bool usingMaxHeFluence;

	/**
	 * The maximum Helium fluence value?
	 */
	double maxHeFluence;

	/**
	 * The amplitude of the flux
	 */
	double heFlux;

public:

	FluxHandler() :
			stepSize(0.0e-16), heFluence(0.0e-16), usingMaxHeFluence(false), maxHeFluence(
					0.0e-16), heFlux(1.0) {
	}

	~FluxHandler() {
	}

	/**
	 * This operation returns the incident flux for a specific cluster composition,
	 * position, and time.
	 * @param compositionVec  The composition of the cluster
	 * @param position        The position of the cluster
	 * @param currentTime     The time
	 * @return incidentFlux   The incident flux at the given position and time of the cluster with
	 * the specified composition
	 */
	virtual double getIncidentFlux(std::vector<int> compositionVec,
			std::vector<double> position, double currentTime);

	/**
	 * Given a specific concentration, position, and time, this operation sets the outgoing
	 * flux to the specified amount.
	 * @param composition  The composition of the cluster
	 * @param position     The position of the cluster
	 * @param time         The time
	 * @return outgoingFlux  The outgoing flux at the given position and time of the cluster with
	 * the specified composition
	 */
	virtual void setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux);

	/**
	 * This operation increments the Helium fluence at the current time step.
	 * @param dt			The length of the time step
	 * @param step			The grid step size
	 * @return				The value of the Helium fluence at the current time step
	 */
	virtual double incrementHeFluence(double dt, double step);

	/**
	 * This operation returns the Helium fluence
	 * @return	The Helium fluence at current time step
	 */
	virtual double getHeFluence() const;

	/**
	 * This operation sets the maximum value of the Helium fluence.
	 * @param fluence	The maximim Helium fluence value
	 */
	virtual void setMaxHeFluence(double fluence);

	/**
	 * This function returns the maximum value of the Helium fluence.
	 */
	virtual double getMaxHeFluence() const;

	/**
	 * This operation sets whether or not the maximum Helium fluence will be used
	 * @param use	If the program should or should not use the max He fluence
	 * @return	True if program will use the max He fluence, and false if it won't
	 */
	virtual bool useMaximumHeFluence();

	/**
	 * This operation sets the factor to change the Helium flux.
	 * @param flux	Helium flux value
	 */
	virtual void setHeFlux(double flux);

	/**
	 * This operation gets the factor that changes the Helium flux.
	 * @return	Helium flux value
	 */
	virtual double getHeFlux() const;

};
//end class FluxHandler

}

#endif
