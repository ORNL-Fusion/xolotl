#ifndef IREACTIONNETWORK_H
#define IREACTIONNETWORK_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include "IReactant.h"

namespace xolotlCore {

class IReactant;
class ProductionReaction;
class DissociationReaction;

/**
 *  This class manages the set of reactants and compound reactants
 *  (combinations of normal reactants). It also manages a set of properties
 *  that describe both.
 */
class IReactionNetwork {

public:
	/**
	 * Nice name for vector of Reactants.
	 */
	using ReactantVector = std::vector<std::reference_wrapper<IReactant> >;

	/**
	 * Nice name for map of reactants, keyed by their composition.
	 */
	using ReactantMap = std::unordered_map<IReactant::Composition, std::unique_ptr<IReactant> >;

	/**
	 * The destructor.
	 */
	virtual ~IReactionNetwork() {
	}

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants.
	 *
	 * @param temp The new temperature
	 */
	virtual void setTemperature(double temp) = 0;

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 *
	 * @return The temperature
	 */
	virtual double getTemperature() const = 0;

	/**
	 * Retrieve the single-species reactant with the given type and size if it
	 * exists in the network or null if not.
	 * Convenience function for get() that takes a
	 * reactant type and composition.
	 *
	 * @param species The reactant's single species.
	 * @param size The size of the reactant.
	 * @return A pointer to the reactant, or nullptr if it does not 
	 * exist in the network.
	 */
	virtual IReactant * get(Species species,
			IReactant::SizeType size) const = 0;

	/**
	 * Retrieve the reactant with the given type and composition if
	 * exists in the network.
	 *
	 * @param type The type of the reactant
	 * @param comp The composition of the reactant
	 * @return A pointer to the reactant of type 'type' and with composition
	 * 'comp.' nullptr if no such reactant exists.
	 */
	virtual IReactant * get(ReactantType type,
			const IReactant::Composition& comp) const = 0;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 *
	 * @return The list of all of the reactants in the network
	 */
	virtual const IReactant::RefVector& getAll() const = 0;

	/**
	 * This operation returns all reactants in the network with the given type.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 *
	 * @param type The reactant or compound reactant type
	 * @return The list of all of the reactants in the network.
	 */
	virtual const ReactantMap& getAll(ReactantType type) const = 0;

	/**
	 * Give the reactant to the network.
	 *
	 * @param reactant The reactant that should be added to the network
	 */
	virtual void add(std::unique_ptr<IReactant> reactant) = 0;

	/**
	 * This operation removes a group of reactants from the network.
	 *
	 * @param reactants The reactants that should be removed.
	 */
	virtual void removeReactants(const ReactantVector& reactants) = 0;

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	virtual void reinitializeNetwork() = 0;

	/**
	 * This method redefines the connectivities for each cluster in the
	 * allReactans vector.
	 */
	virtual void reinitializeConnectivities() = 0;

	/**
	 * This operation returns the size or number of reactants in the network.
	 *
	 * @return The number of reactants in the network
	 */
	virtual int size() const = 0;

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() const = 0;

	/**
	 * This operation returns the list (vector) of each reactant in the network.
	 *
	 * @return The list of compositions
	 */
	virtual std::vector< std::vector <int> > getCompositionList() const = 0;

	/**
	 * Add a production reaction to the network.
	 *
	 * @param reaction The reaction that should be added to the network
	 * @return The reaction that is now in the network
	 */
	virtual ProductionReaction& add(
			std::unique_ptr<ProductionReaction> reaction) = 0;

	/**
	 * Add a dissociation reaction to the network.
	 *
	 * @param reaction The reaction that should be added to the network
	 * @return The reaction that is now in the network
	 */
	virtual DissociationReaction& add(
			std::unique_ptr<DissociationReaction> reaction) = 0;

	/**
	 * This operation fills an array of doubles with the concentrations of all
	 * of the reactants in the network.
	 *
	 * @param concentrations The array that will be filled with the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. If the array is too small to hold the concentrations, SIGSEGV will
	 * be thrown.
	 */
	virtual void fillConcentrationsArray(double * concentrations) = 0;

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	virtual void updateConcentrationsFromArray(double * concentrations) = 0;

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * @param diagFill The pointer to the vector where the connectivity information is kept
	 */
	virtual void getDiagonalFill(int *diagFill) = 0;

	/**
	 * Get the total concentration of atoms in the network.
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	virtual double getTotalAtomConcentration(int i = 0) = 0;

	/**
	 * Get the total concentration of atoms contained in bubbles in the network.
	 *
	 * @param i Index to switch between the different types of atoms
	 * @return The total concentration
	 */
	virtual double getTotalTrappedAtomConcentration(int i = 0) = 0;

	/**
	 * Get the total concentration of vacancies in the network.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalVConcentration() = 0;

	/**
	 * Get the total concentration of material interstitials in the network.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalIConcentration() = 0;

	/**
	 * Calculate all the rate constants for the reactions and dissociations of the network.
	 * Need to be called only when the temperature changes.
	 */
	virtual void computeRateConstants() = 0;

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 */
	virtual void computeAllFluxes(double *updatedConcOffset) = 0;

    /**
     * Determine the number of partials for each cluster
     * and their starting locations within the vectors used
     * when computing partials.
     *
     * @param size Number of partials for each cluster.
     * @param startingIdx Starting index of items owned by each reactant
     *      within the partials values array and the indices array.
     * @return Total number of partials for all clusters.
     */
    virtual size_t initPartialsSizes(std::vector<int>& size,
                                    std::vector<size_t>& startingIdx) const = 0;

    /**
     * Initialize the indexing used for computing partial derivatives
     * for all clusters and their momenta.
     *
     * @param size Number of partials for each cluster.
     * @param startingIdx Starting index of items owned by each reactant
     *      within the partials values array and the indices array.
     * @param indices The indices of the clusters for the partial derivatives.
     */
    virtual void initPartialsIndices(const std::vector<int>& size,
                                    const std::vector<size_t>& startingIdx,
                                    std::vector<int>& indices) const = 0;

    /**
     * Compute the partial derivatives generated by all the reactions
     * for all the clusters and their momentum.
     *
     * @param size Number of partials for each cluster.
     * @param startingIdx Starting index of items owned by each reactant
     *      within the partials values array and the indices array.
     * @param indices The indices of the clusters for the partial derivatives.
     * @param vals The values of partials for the reactions
     */
    virtual void computeAllPartials(const std::vector<int>& size,
                                    const std::vector<size_t>& startingIdx,
                                    const std::vector<int>& indices,
                                    std::vector<double>& vals) const = 0;

	/**
	 * This operation returns the biggest production rate in the network.
	 *
	 * @return The biggest rate
	 */
	virtual double getBiggestRate() const = 0;

	/**
	 * Are dissociations enabled?
	 *
	 * @returns true if reactions are enabled, false otherwise.
	 */
	virtual bool getDissociationsEnabled() const = 0;

	/**
	 * Find maximum cluster size currently in the network
	 * for the given reactant type.
	 *
	 * @param rtype Reactant type of interest.
	 * @return Maximum size of cluster of type rtype currently in network.
	 */
	virtual IReactant::SizeType getMaxClusterSize(ReactantType rtype) const = 0;

	/**
	 * Dump a representation of the network to the given output stream.
	 *
	 * @param os Output stream on which to write network description.
	 */
	virtual void dumpTo(std::ostream& os) const = 0;
};

}

#endif
