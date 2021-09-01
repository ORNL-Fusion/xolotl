#ifndef ALPHAZRFITFLUXHANDLER_H
#define ALPHAZRFITFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/ZrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident
 * fluxes for an alpha Zr material.
 */
class AlphaZrFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Not actually used
		return 0.0;
	}

public:
	/**
	 * The constructor
	 */
	AlphaZrFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~AlphaZrFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Only defined in 0D
		if (xGrid.size() == 0) {
			// Add an empty vector
			std::vector<double> tempVector;
			incidentFluxVec.push_back(tempVector);
		}

		using NetworkType = network::ZrReactionNetwork;
		auto zrNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single interstitial cluster
		// here
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::I] = 1;
		auto cluster = zrNetwork->findCluster(comp, plsm::onHost);
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe single interstitial cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
        
        
        comp[NetworkType::Species::I] = 2;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 2 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 3;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 3 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 4;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 4 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 5;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 5 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 6;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 6 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 7;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 7 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 8;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 8 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::I] = 9;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 9 interstitial cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        

		// Look for vacancy now
		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::V] = 1;
		cluster = zrNetwork->findCluster(comp, plsm::onHost);
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
        
        
        comp[NetworkType::Species::V] = 2;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error(
                                     "\nThe 2 vacancy cluster is not present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::V] = 3;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 3 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        
        comp[NetworkType::Species::V] = 4;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 4 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::V] = 5;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 5 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::V] = 6;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 6 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::V] = 7;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 7 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::V] = 8;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 8 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
        comp[NetworkType::Species::V] = 9;
        cluster = zrNetwork->findCluster(comp, plsm::onHost);
        if (cluster.getId() == NetworkType::invalidIndex()) {
            throw std::runtime_error("\nThe 9 vacancy cluster is not "
                                     "present in the network, "
                                     "cannot use the flux option!");
        }
        fluxIndices.push_back(cluster.getId());
        
		return;
	}
         
	/**
	 * \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		// Define only for a 0D case
		if (incidentFluxVec[0].size() == 0) {
   
			updatedConcOffset[fluxIndices[0]] += 7.15306e-07; // I1 (7.153e+14 cm^-3 s^-1)
			updatedConcOffset[fluxIndices[1]] += 4.37883e-08; // I2
            updatedConcOffset[fluxIndices[2]] += 2.43045e-08; // I3
            updatedConcOffset[fluxIndices[3]] += 1.61015e-08; // I4
            updatedConcOffset[fluxIndices[4]] += 4.66543e-09; // I5
            updatedConcOffset[fluxIndices[5]] += 5.13252e-09; // I6
            updatedConcOffset[fluxIndices[6]] += 5.13252e-09; // I7
            updatedConcOffset[fluxIndices[7]] += 5.13252e-09; // I8
            updatedConcOffset[fluxIndices[8]] += 5.13252e-09; // I9
            updatedConcOffset[fluxIndices[9]] += 5.18329e-07; // V1
            updatedConcOffset[fluxIndices[10]] += 5.32574e-08; // V2
            updatedConcOffset[fluxIndices[11]] += 1.28488e-08; // V3
            updatedConcOffset[fluxIndices[12]] += 2.01684e-08; // V4
            updatedConcOffset[fluxIndices[13]] += 1.51447e-08; // V5
            updatedConcOffset[fluxIndices[14]] += 2.91043e-09; // V6
            updatedConcOffset[fluxIndices[15]] += 2.91043e-09; // V7
            updatedConcOffset[fluxIndices[16]] += 2.91043e-09; // V8
            updatedConcOffset[fluxIndices[17]] += 2.91043e-09; // V9
            
		}

		else {
			throw std::runtime_error(
				"\nThe alpha Zr problem is not defined for more than 0D!");
		}

		return;
	}
};
// end class AlphaZrFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
