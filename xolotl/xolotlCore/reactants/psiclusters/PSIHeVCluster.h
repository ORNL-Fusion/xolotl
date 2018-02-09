#ifndef PSIHEVCLUSTER_H
#define PSIHEVCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  A cluster composed of helium and vacancies
 */
class PSIHeVCluster: public PSICluster {

private:
	// TODO do we need to keep these species counts here,
	// since they are in the composition?

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of atomic vacancies in this cluster.
	int numV;

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
	PSIHeVCluster() = delete;

	/**
	 * The constructor. All PSIHeVClusters must be initialized with a map
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
	PSIHeVCluster(int numHe, int numV, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(numHe, numV)), numHe(
					numHe), numV(numV) {
		// Set the cluster size as the sum of
		// the number of Helium and Vacancies
		size = numHe + numV;

		// Update the composition map
		composition[toCompIdx(Species::He)] = numHe;
		composition[toCompIdx(Species::V)] = numV;

		// Set the typename appropriately
		type = ReactantType::HeV;

		// Compute the reaction radius
		reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
				+ pow(
						(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0)
								* numV) / (8.0 * xolotlCore::pi), (1.0 / 3.0))
				- pow(
						(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
								/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

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
	PSIHeVCluster(PSIHeVCluster &other) = delete;

	//! Destructor
	~PSIHeVCluster() {
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
//end class PSIHeVCluster

} /* end namespace xolotlCore */
#endif
