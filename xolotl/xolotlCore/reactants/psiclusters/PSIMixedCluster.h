#ifndef PSIMIXEDCLUSTER_H
#define PSIMIXEDCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  A cluster composed of impurities and vacancies
 */
class PSIMixedCluster: public PSICluster {

private:
	// TODO do we need to keep these species counts here,
	// since they are in the composition?

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of deuterium atoms in this cluster.
	int numD;

	//! The number of tritium atoms in this cluster.
	int numT;

	//! The number of atomic vacancies in this cluster.
	int numV;

	static std::string buildName(IReactant::SizeType nHe,
			IReactant::SizeType nD, IReactant::SizeType nT,
			IReactant::SizeType nV) {
		std::stringstream nameStream;
		nameStream << "He_" << nHe << "D_" << nD << "T_" << nT << "V_" << nV;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSIMixedCluster() = delete;

	/**
	 * The constructor. All PSIMixedClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {He,V}.
	 *
	 * @param numHe The number of helium atoms in this cluster
	 * @param numD The number of deuterium in this cluster
	 * @param numT The number of tritium in this cluster
	 * @param numV The number of vacancies in this cluster
	 * @param _network The network the cluster will belong to.
	 * @param registry The performance handler registry
	 */
	PSIMixedCluster(int numHe, int numD, int numT, int numV,
			IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(numHe, numD, numT, numV)), numHe(
					numHe), numD(numD), numT(numT), numV(numV) {
		// Set the cluster size as the sum of
		// the number of Helium, Hydrogen, Vacancies
		size = numHe + numD + numT + numV;

		// Update the composition map
		composition[toCompIdx(Species::He)] = numHe;
		composition[toCompIdx(Species::D)] = numD;
		composition[toCompIdx(Species::T)] = numT;
		composition[toCompIdx(Species::V)] = numV;

		// Set the typename appropriately
		type = ReactantType::PSIMixed;

		// Compute the reaction radius
		double latticeParam = network.getLatticeParameter();
		reactionRadius = (sqrt(3.0) / 4.0) * latticeParam
				+ pow(
						(3.0 * pow(latticeParam, 3.0) * numV)
								/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
				- pow((3.0 * pow(latticeParam, 3.0)) / (8.0 * xolotlCore::pi),
						(1.0 / 3.0));

		// Bounds on He, D, T, and V
		bounds[0] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(numHe),
				static_cast<IReactant::SizeType>(numHe + 1));
		bounds[1] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(numD),
				static_cast<IReactant::SizeType>(numD + 1));
		bounds[2] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(numT),
				static_cast<IReactant::SizeType>(numT + 1));
		bounds[3] = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(numV),
				static_cast<IReactant::SizeType>(numV + 1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	PSIMixedCluster(PSIMixedCluster &other) = delete;

	//! Destructor
	~PSIMixedCluster() {
	}

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const override {
		return true;
	}

	/**
	 * Add grid points to the vector of diffusion coefficients or remove
	 * them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	void addGridPoints(int i) override {
		// Don't do anything
		return;
	}

	/**
	 * This operation sets the temperature at which the reactant currently
	 * exists. Temperature-dependent quantities are recomputed when this
	 * operation is called, so the temperature should always be set first.
	 *
	 * @param temp The new cluster temperature
	 * @param i The location on the grid
	 */
	void setTemperature(double temp, int i) override {
		// Don't do anything
		return;
	}

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @param i The position on the grid
	 * @return The diffusion coefficient
	 */
	double getDiffusionCoefficient(int i) const override {
		return 0.0;
	}

};
//end class PSIMixedCluster

} /* end namespace xolotlCore */
#endif
