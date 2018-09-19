#ifndef FEHEVCLUSTER_H
#define FEHEVCLUSTER_H

// Includes
#include "FeCluster.h"
#include <string>
#include <map>
#include <Constants.h>

namespace xolotlCore {

/**
 *  A cluster composed of helium and vacancies
 */
class FeHeVCluster: public FeCluster {

private:

	static std::string buildName(IReactant::SizeType nHe,
			IReactant::SizeType nV) {
		std::stringstream nameStream;
		nameStream << "He_" << nHe << "V_" << nV;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FeHeVCluster() = delete;

	/**
	 * The constructor. All HeVClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {He,V}.
	 *
	 * @param numHe The number of helium atoms in this cluster
	 * @param numV The number of vacancies in this cluster
	 * @param _network The network the cluster will belong to.
	 * @param registry The performance handler registry
	 */
	FeHeVCluster(int numHe, int numV, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			FeCluster(_network, registry, buildName(numHe, numV)) {
		// Set the cluster size as the sum of
		// the number of Helium and Vacancies
		size = numHe + numV;

		// Update the composition map
		composition[toCompIdx(Species::He)] = numHe;
		composition[toCompIdx(Species::V)] = numV;

		// Set the typename appropriately
		type = ReactantType::HeV;

		// Compute the reaction radius
		// It is the same formula for HeV clusters
		reactionRadius = xolotlCore::ironLatticeConstant
				* pow((3.0 * numV) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(numHe),
				static_cast<IReactant::SizeType>(numHe + 1));
		vBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(numV),
				static_cast<IReactant::SizeType>(numV + 1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeHeVCluster(FeHeVCluster &other) = delete;

	//! Destructor
	~FeHeVCluster() {
	}

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {
		return true;
	}

};
//end class FeHeVCluster

} /* end namespace xolotlCore */
#endif
