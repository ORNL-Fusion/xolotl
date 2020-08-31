#ifndef PULSEDFITFLUXHANDLER_H
#define PULSEDFITFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident V
 * and I flux for tungsten material with a flux amplitude being intermittent.
 */
class PulsedFitFluxHandler : public FluxHandler
{
private:
	/**
	 * Total time length of the pulse
	 */
	double deltaTime;

	/**
	 * Proportion of the total time where the flux amplitude is not 0.0
	 */
	double alpha;

	/**
	 * Parameters for the gaussian profile
	 */
	double mu = 2000.0;
	double sigma = 100.0;

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double
	FitFunction(double x)
	{
		// Compute the polynomial fit
		double value = exp(-pow((x - mu) / (sqrt(2.0) * sigma), 2.0));

		return std::max(value, 0.0);
	}

public:
	/**
	 * The constructor
	 */
	PulsedFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options),
		deltaTime(0.0),
		alpha(0.0)
	{
	}

	/**
	 * The Destructor
	 */
	~PulsedFitFluxHandler()
	{
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (util::equal(fluxAmplitude, 0.0) && !useTimeProfile) {
			return;
		}

		// Set the flux index corresponding the the single vacancy cluster here
		using NetworkType = network::IPSIReactionNetwork;
		using AmountType = NetworkType::AmountType;

		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto numSpecies = psiNetwork->getSpeciesListSize();
		auto specIdV = psiNetwork->getVacancySpeciesId();
		auto specIdI = psiNetwork->getInterstitialSpeciesId();

		// Set the flux index corresponding the the single helium cluster here
		auto comp = std::vector<AmountType>(numSpecies, 0);
		comp[specIdV()] = 1;
		auto clusterId = psiNetwork->findClusterId(comp);
		// Check that the vacancy cluster is present in the network
		if (clusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(clusterId);

		// Set the flux index corresponding the the single interstitial cluster
		// here Initialize the composition
		comp.assign(numSpecies, 0);
		comp[specIdI()] = 1;
		clusterId = psiNetwork->findClusterId(comp);
		// Check that the interstitial cluster is present in the network
		if (clusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe single interstitial cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(clusterId);
	}

	/**
	 * This operation computes the flux due to incoming particles at a given
	 * grid point. \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		// Check in which phase of the pulse we are
		int cycle = currentTime / deltaTime;
		// The flux is 0.0 after alpha * deltaTime
		if (currentTime - ((double)cycle * deltaTime) > alpha * deltaTime ||
			util::equal(deltaTime, 0.0) || util::equal(alpha, 0.0))
			return;

		// Update the concentration array
		updatedConcOffset[fluxIndices[0]] +=
			incidentFluxVec[0][xi - surfacePos]; // V
		updatedConcOffset[fluxIndices[1]] +=
			incidentFluxVec[0][xi - surfacePos]; // I

		return;
	}

	/**
	 * This operation sets the time of the pulse.
	 * \see IFluxHandler.h
	 */
	void
	setPulseTime(double time)
	{
		deltaTime = time;
		return;
	}

	/**
	 * This operation sets proportion of the pulse that is on.
	 * \see IFluxHandler.h
	 */
	void
	setProportion(double a)
	{
		alpha = a;
		return;
	}
};
// end class PulsedFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
