#ifndef FLUXHANDLER_H
#define FLUXHANDLER_H

#include "IFluxHandler.h"
#include <vector>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * flux calculations.
 */
class FluxHandler: public IFluxHandler {

protected:

	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position).
	 */
	std::vector<double> incidentFluxVec;

	/**
	 * Step size between each grid point in the x direction.
	 */
	double stepSize;

	/**
	 * Helium fluence.
	 */
	double heFluence;

	/**
	 * The amplitude of the flux.
	 */
	double heFlux;

	/**
	 * Are we using a time profile for the amplitude of the helium incoming flux?
	 */
	bool useTimeProfile;

	/**
	 * Value of the fit function integrated on the grid.
	 */
	double normFactor;

	/**
	 * Vector to hold the time read from the input
	 * time profile file.
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the amplitude read from the input
	 * time profile file.
	 */
	std::vector<double> amplitude;

	/**
	 * Function that calculates the flux at a given position x (in nm).
	 * It needs to be implemented by the daughter classes.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	virtual double FitFunction(double x) {return 0.0;}

	/**
	 * This method returns the value of the helium incident flux amplitude at the
	 * given time when a time profile is used.
	 *
	 * @param currentTime The time
	 * @return The value of the helium flux at this time
	 */
	double getAmplitude(double currentTime) const;

	/**
	 * This method recomputes the values of the incident flux vector when
	 * a time profile is given.
	 */
	void recomputeFluxHandler();

public:

	FluxHandler();

	~FluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
     * \see IFluxHandler.h
	 */
	virtual void initializeFluxHandler(int nx, double hx);

	/**
	 * This method reads the values on the time profile file and store them in the
	 * time and amplitude vectors.
     * \see IFluxHandler.h
	 */
	void initializeTimeProfile(const std::string& fileName);

	/**
	 * This operation returns the incident flux vector.
     * \see IFluxHandler.h
	 */
	virtual std::vector<double> getIncidentFluxVec(double currentTime);

	/**
	 * This operation increments the helium fluence at the current time step.
     * \see IFluxHandler.h
	 */
	virtual void incrementHeFluence(double dt);

	/**
	 * This operation returns the helium fluence.
     * \see IFluxHandler.h
	 */
	virtual double getHeFluence() const;

	/**
	 * This operation sets the factor to change the intensity of the helium flux.
     * \see IFluxHandler.h
	 */
	virtual void setHeFlux(double flux);

	/**
	 * This operation gets the factor that changes the helium flux intensity.
     * \see IFluxHandler.h
	 */
	virtual double getHeFlux() const;

};
//end class FluxHandler

}

#endif
