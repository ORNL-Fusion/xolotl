#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <IReactionNetwork.h>
#include <memory>

namespace xolotlCore {

namespace HDF5Utils {

	/**
	 * Create the HDF5 file with the needed structure.
	 *
	 * @param fileName The name of the file to create
	 */
	void initializeFile(const std::string& fileName);

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
	 * @param fileName The name of the file from which the network will be taken
	 */
	void fillNetwork(const std::string& fileName);

	/**
	 * Add a concentration subgroup for the given time step to the HDF5 file.
	 *
	 * @param timeStep The number of the time step
	 * @param time The physical time at this time step
	 * @param previousTime The physical time at the previous time step
	 * @param deltaTime The physical length of the time step
	 */
	void addConcentrationSubGroup(int timeStep, double time,
			double previousTime, double deltaTime);

	/**
	 * Write the surface position as an attribute of the
	 * concentration subgroup.
	 *
	 * @param timeStep The number of the time step
	 * @param iSurface The index of the surface position
	 * @param nInter The quantity of interstitial at each surface position
	 * @param previousFlux The previous I flux at each surface position
	 */
	void writeSurface1D(int timeStep, int iSurface,
			double nInter, double previousFlux);

	/**
	 * Write the surface positions as a dataset of the
	 * concentration subgroup.
	 *
	 * @param timeStep The number of the time step
	 * @param iSurface The indices of the surface position
	 * @param nInter The quantity of interstitial at each surface position
	 * @param previousFlux The previous I flux at each surface position
	 */
	void writeSurface2D(int timeStep, std::vector<int> iSurface,
			std::vector<double> nInter, std::vector<double> previousFlux);

	/**
	 * Write the surface positions as a dataset of the
	 * concentration subgroup.
	 *
	 * @param timeStep The number of the time step
	 * @param iSurface The indices of the surface position
	 * @param nInter The quantity of interstitial at each surface position
	 * @param previousFlux The previous I flux at each surface position
	 */
	void writeSurface3D(int timeStep,
			std::vector< std::vector<int> > iSurface,
			std::vector< std::vector<double> > nInter,
			std::vector< std::vector<double> > previousFlux);

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
	 * @param concArray The array of concentration at a grid point
	 * @param i The index of the position on the grid on the x direction
	 * @param j The index of the position on the grid on the y direction
	 * @param k The index of the position on the grid on the z direction
	 */
	void fillConcentrations(double concArray[][2],
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
	 * Read the previous time from the concentration group of a HDF5 file.
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The physical time at the previous timestep
	 */
	double readPreviousTime(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the surface position from the concentration group of a HDF5 file in
	 * the case of a 1D grid (one surface position).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The index of the surface position
	 */
	int readSurface1D(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the surface position from the concentration group of a HDF5 file in
	 * the case of a 2D grid (a vector of surface positions).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The vector of indices of the surface position
	 */
	std::vector<int> readSurface2D(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the surface position from the concentration group of a HDF5 file in
	 * the case of a 3D grid (a vector of vector of surface positions).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The vector of vector of indices of the surface position
	 */
	std::vector< std::vector<int> > readSurface3D(const std::string& fileName,
			int lastTimeStep);

	/**
	 * Read the quantity of interstitial at the surface from the concentration
	 * group of a HDF5 file in the case of a 1D grid (one surface position).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The number of interstitial
	 */
	double readNInterstitial1D(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the quantity of interstitial at each surface position
	 * from the concentration group of a HDF5 file in
	 * the case of a 2D grid (a vector of surface positions).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The vector of interstitial quantity of the surface position
	 */
	std::vector<double> readNInterstitial2D(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the quantity of interstitial at each surface position
	 * from the concentration group of a HDF5 file in
	 * the case of a 3D grid (a vector of vector of surface positions).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The vector of interstitial quantity of the surface position
	 */
	std::vector< std::vector<double> > readNInterstitial3D(const std::string& fileName,
			int lastTimeStep);

	/**
	 * Read the previous interstitial flux at the surface from the concentration
	 * group of a HDF5 file in the case of a 1D grid (one surface position).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The number of interstitial
	 */
	double readPreviousIFlux1D(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the previous interstitial flux at each surface position
	 *  from the concentration group of a HDF5 file in
	 * the case of a 2D grid (a vector of surface positions).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The vector of previous flux at each surface position
	 */
	std::vector<double> readPreviousIFlux2D(const std::string& fileName, int lastTimeStep);

	/**
	 * Read the previous interstitial flux at each surface position
	 * from the concentration group of a HDF5 file in
	 * the case of a 3D grid (a vector of vector of surface positions).
	 *
	 * @param fileName The name of the file to read from
	 * @param lastTimeStep The value of the last written time step
	 * @return The vector of previous flux at each surface position
	 */
	std::vector< std::vector<double> > readPreviousIFlux3D(const std::string& fileName,
			int lastTimeStep);

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

}

} /* namespace xolotlCore */
#endif
