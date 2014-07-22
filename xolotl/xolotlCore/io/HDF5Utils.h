#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

namespace HDF5Utils {

	/**
	 * Create the HDF5 file with the needed structure.
	 * @param fileName The name of the file to create.
	 * @param networkSize The total number of cluster in the network.
	 * @param gridSize The total number of grid points.
	 */
	void initializeFile(std::string fileName, int networkSize, int gridSize);

	/**
	 * Open the already existing HDF5 file.
	 * @param fileName The name of the file to open.
	 */
	void openFile(std::string fileName);

	/**
	 * Fill the header.
	 * @param physicalDim The physical length of the material on which one is solving the ADR equation.
	 * @param refinement The refinement of the grid.
	 */
	void fillHeader(int physicalDim, int refinement);

	/**
	 * Fill the network dataset.
	 * @param network The network of clusters.
	 */
	void fillNetwork(std::shared_ptr<PSIClusterReactionNetwork> network);

	/**
	 * Add a concentration subgroup for the given time step to the HDF5 file.
	 * @param timeStep The number of the time step.
	 * @param networkSize The total number of cluster in the network.
	 * @param time The physical time at this time step.
	 * @param deltaTime The physical length of the time step.
	 */
	void addConcentrationSubGroup(int timeStep, int networkSize, double time,
			double deltaTime);

	/**
	 * Fill the concentration dataset at a specific grid point.
	 * @param index The index of the position on the grid.
	 * @param size The size of the dataset to create.
	 */
	void addConcentrationDataset(int index, int size);

	/**
	 * Fill the concentration dataset at a specific grid point.
	 * @param concVector The vector of concentration at a grid point.
	 * @param index The index of the position on the grid.
	 */
	void fillConcentrations(std::vector< std::vector<double> > concVector, int index);

	/**
	 * Close the file for the first time after creating it.
	 */
	void finalizeFile();

	/**
	 * Close the file when it had been opened by openFile().
	 */
	void closeFile();

	/**
	 * Read the header of a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @param physicalDim The physical length of the material to be changed.
	 */
	void readHeader(std::string fileName, int & physicalDim);

	/**
	 * Check if the file contains a valid concentration group.
	 * @param fileName The name of the file to read from.
	 * @param lastTimeStep The value of the last written time step to be changed.
	 * @return True is the file contains a valid concentration group.
	 */
	bool hasConcentrationGroup(std::string fileName, int & lastTimeStep);

	/**
	 * Read the times from the concentration group of a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @param lastTimeStep The value of the last written time step.
	 * @param time The physical time to be changed.
	 * @param deltaTime The time step length to be changed.
	 */
	void readTimes(std::string fileName, int lastTimeStep, double & time, double & deltaTime);

	/**
	 * Read the network from a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @return The vector of vector which contain the network dataset.
	 */
	std::vector< std::vector <double> > readNetwork(std::string fileName);

	/**
	 * Read the i-th grid point concentrations from a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @param lastTimeStep The value of the last written time step.
	 * @param i The index of the grid point.
	 * @return The vector of concentrations.
	 */
	std::vector< std::vector<double> > readGridPoint(std::string fileName,
			int lastTimeStep, int i);

};

} /* namespace xolotlCore */
#endif
