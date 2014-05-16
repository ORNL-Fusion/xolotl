#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

namespace HDF5Utils {

	/**
	 * Create the HDF5 file with the needed structure.
	 * @param timeStep The number of the time step.
	 * @param networkSize The total number of cluster in the network.
	 */
	void initializeFile(int timeStep, int networkSize);

	/**
	 * Fill the header.
	 * @param physicalDim The physical lenght of the material on which one is solving the ADR equation.
	 * @param refinement The refinement of the grid.
	 * @param time The physical time at this time step.
	 * @param deltaTime The physical length of the time step.
	 */
	void fillHeader(int physicalDim, int refinement, double time, double deltaTime);

	/**
	 * Fill the network dataset.
	 * @param network The network of clusters.
	 */
	void fillNetwork(std::shared_ptr<PSIClusterReactionNetwork> network);

	/**
	 * Fill the concentration dataset at a specific grid point.
	 * @param concArray The vector of concentration at a grid point.
	 * @param position The physical position on the grid.
	 */
	void fillConcentrations(double * concArray, int index, double position);

	/**
	 * Add the data to the file and close it.
	 */
	void finalizeFile();

	/**
	 * Read the network from a HDF5 file.
	 * @param fileName The name of the file to read from.
	 */
	void readNetwork(std::string fileName);

};

} /* namespace xolotlCore */
#endif
