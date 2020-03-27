#ifndef UZRXEVCLUSTER_H
#define UZRXEVCLUSTER_H

// Includes
#include "UZrCluster.h"
#include <string>
#include <map>
#include <Constants.h>

namespace xolotlCore {

/**
 *  A cluster composed of xenon and vacancies
 */
class UZrXeVCluster: public UZrCluster {

private:
	//! The number of xenon atoms in this cluster.
	int numXe;

	//! The number of atomic vacancies in this cluster.
	int numV;

	static std::string buildName(IReactant::SizeType nXe,
			IReactant::SizeType nV) {
		std::stringstream nameStream;
		nameStream << "Xe_" << nXe << "V_" << nV;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	UZrXeVCluster() = delete;

	/**
	 * The constructor. All XeVClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {Xe,V}.
	 *
	 * @param numXe The number of xenon atoms in this cluster
	 * @param numV The number of vacancies in this cluster
	 * @param _network The network the cluster will belong to.
	 * @param registry The performance handler registry
	 */
	UZrXeVCluster(int numXe, int numV, IReactionNetwork &_network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			UZrCluster(_network, registry, buildName(numXe, numV)), numXe(
					numXe), numV(numV) {
		// Set the cluster size as the sum of
		// the number of Xenon and Vacancies
		size = numXe + numV;

		// Update the composition map
		composition[toCompIdx(Species::Xe)] = numXe;
		composition[toCompIdx(Species::V)] = numV;

		// Set the typename appropriately
		type = ReactantType::XeV;

		// TODO: Compute the reaction radius
		// You can look at FeHeVCluster or PSIMixedCluster
		reactionRadius = network.getLatticeParameter();
		;

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	UZrXeVCluster(UZrXeVCluster &other) = delete;

	//! Destructor
	~UZrXeVCluster() {
	}

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * Xe and V.
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
//end class UZrXeVCluster

} /* end namespace xolotlCore */
#endif
