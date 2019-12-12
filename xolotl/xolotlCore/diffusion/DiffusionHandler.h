#ifndef DIFFUSIONHANDLER_H
#define DIFFUSIONHANDLER_H

// Includes
#include "IDiffusionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile clusters. It needs to have
 * subclasses implementing the compute diffusion methods.
 */
class DiffusionHandler: public IDiffusionHandler {
protected:

	//! Collection of diffusing clusters.
	IReactant::ConstRefVector diffusingClusters;

	//! If we want the hydrogen desorption
	bool isDesorbing;

public:

	//! The Constructor
	DiffusionHandler() :
			isDesorbing(false) {
	}

	//! The Destructor
	~DiffusionHandler() {
	}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped it
	 * won't be possible to set the partial derivatives for the diffusion.
	 *
	 * The value 1 is set in ofillMap if a cluster has a non zero diffusion coefficient.
	 *
	 * @param network The network
	 * @param ofillMap Map of connectivity for diffusing clusters.
	 */
	virtual void initializeOFill(const IReactionNetwork& network,
			IReactionNetwork::SparseFillMap& ofillMap) override {

		int dof = network.getDOF();

		// Clear the index vector
		diffusingClusters.clear();

		// Consider each cluster.
		for (IReactant const& currReactant : network.getAll()) {

			auto const& cluster = static_cast<IReactant const&>(currReactant);

			// Get its diffusion coefficient
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (xolotlCore::equal(diffFactor, 0.0))
				continue;

			// Note that cluster is diffusing.
			diffusingClusters.emplace_back(cluster);

			// Get its id
			int index = cluster.getId() - 1;
			// Set the ofill value to 1 for this cluster
			ofillMap[index].emplace_back(index);
		}

		return;
	}

	/**
	 * Get the total number of diffusing clusters in the network.
	 *
	 * @return The number of diffusing clusters
	 */
	int getNumberOfDiffusing() const override {
		return diffusingClusters.size();
	}

	/**
	 * Set if we want the hydrogen desorption.
	 *
	 * @param desorb True if we want to desorb
	 */
	void setDesorption(bool desorb) override {
		isDesorbing = desorb;
	}

};
//end class DiffusionHandler

} /* end namespace xolotlCore */
#endif
