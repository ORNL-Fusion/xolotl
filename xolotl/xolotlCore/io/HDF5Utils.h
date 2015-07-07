#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

namespace HDF5Utils {

	/**
	 * Create the HDF5 file with the needed structure.
	 *
	 * @param fileName The name of the file to create
	 * @param networkSize The total number of cluster in the network
	 */
	void initializeFile(const std::string& fileName, int networkSize);

	/**
	 * Open the already existing HDF5 file.
	 *
	 * @param fileName The name of the file to open
	 */
	void openFile(const std::string& fileName);

	/**
	 * Fill the header with the number of points and step size in
	 * each direction.
	 *
	 * @param nx The number of grid points in the x direction (depth)
	 * @param hx The step size in the x direction
	 * @param ny The number of grid points in the y direction
	 * @param hy The step size in the y direction
	 * @param nz The number of grid points in the z direction
	 * @param hz The step size in the z direction
	 */
	void fillHeader(int nx, double hx, int ny = 0,
			double hy = 0.0, int nz = 0, double hz = 0.0);

	/**
	 * Fill the network dataset.
	 *
	 * @param network The network of clusters
	 */
	void fillNetwork(PSIClusterReactionNetwork *network);

	/**
	 * Add a concentration subgroup for the given time step to the HDF5 file.
	 *
	 * @param timeStep The number of the time step
	 * @param networkSize The total number of cluster in the network
	 * @param time The physical time at this time step
	 * @param deltaTime The physical length of the time step
	 */
	void addConcentrationSubGroup(int timeStep, int networkSize, double time,
			double deltaTime);

	/**
	 * Add the concentration dataset at a specific grid point.
	 *
	 * @param size The size of the dataset to create
	 * @param i The index of the position on the grid on the x direction
	 * @param j The index of the position on the grid on the y direction
	 * @param k The index of the position on the grid on the z direction
	 */
	void addConcentrationDataset(int size, int i, int j = -1, int k = -1);

	/**
	 * Fill the concentration dataset at a specific grid point.
	 *
	 * @param concVector The vector of concentration at a grid point
	 * @param i The index of the position on the grid on the x direction
	 * @param j The index of the position on the grid on the y direction
	 * @param k The index of the position on the grid on the z direction
	 */
	void fillConcentrations(const std::vector< std::vector<double> >& concVector,
			int i, int j = -1, int k = -1);

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
	 *
	 * @param fileName The name of the file to read from
	 * @param nx The number of grid points in the x direction (depth)
	 * @param hx The step size in the x direction
	 * @param ny The number of grid points in the y direction
	 * @param hy The step size in the y direction
	 * @param nz The number of grid points in the z direction
	 * @param hz The step size in the z direction
	 */
	void readHeader(const std::string& fileName, int &nx, double &hx, int &ny,
			double &hy, int &nz, double &hz);

	/**
	 * Check if the file contains a valid concentration group.
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step to be changed
	 * @return True if the file contains a valid concentration group
	 */
	bool hasConcentrationGroup(const std::string& fileName, int &lastTimeStep);

	/**
	 * Read the times from the concentration group of a HDF5 file.
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @param time The physical time to be changed
	 * @param deltaTime The time step length to be changed
	 */
	void readTimes(const std::string& fileName, int lastTimeStep, double &time,
			double &deltaTime);

	/**
	 * Read the network from a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @return The vector of vector which contain the network dataset.
	 */
	std::vector< std::vector <double> > readNetwork(const std::string& fileName);

	/**
	 * Read the (i,j,k)-th grid point concentrations from a HDF5 file.
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @param i The index of the grid point on the x axis
	 * @param j The index of the grid point on the y axis
	 * @param k The index of the grid point on the z axis
	 * @return The vector of concentrations
	 */
	std::vector< std::vector<double> > readGridPoint(const std::string& fileName,
			int lastTimeStep, int i, int j = -1, int k = -1);

};

} /* namespace xolotlCore */
#endif
