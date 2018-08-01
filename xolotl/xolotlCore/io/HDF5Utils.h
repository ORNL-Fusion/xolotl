#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <IReactionNetwork.h>
#include <memory>
#include <set>

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
 * @param grid The grid points in the x direction (depth)
 * @param ny The number of grid points in the y direction
 * @param hy The step size in the y direction
 * @param nz The number of grid points in the z direction
 * @param hz The step size in the z direction
 */
void fillHeader(std::vector<double>& grid, int ny = 0, double hy = 0.0, int nz =
		0, double hz = 0.0);

/**
 * Fill the list of cluster compositions.
 *
 * @param compVec The vector of composition
 */
void fillNetworkComp(std::vector<std::vector<int> > compVec);

/**
 * Fill the network dataset.
 *
 * @param fileName The name of the file from which the network will be taken
 * @param network The network to write
 */
void fillNetwork(const std::string& fileName, IReactionNetwork& network);

/**
 * Add a concentration subgroup for the given time step to the HDF5 file.
 *
 * @param timeStep The number of the time step
 * @param time The physical time at this time step
 * @param previousTime The physical time at the previous time step
 * @param deltaTime The physical length of the time step
 */
void addConcentrationSubGroup(int timeStep, double time, double previousTime,
		double deltaTime);

/**
 * Write the surface position as an attribute of the
 * concentration subgroup.
 *
 * @param timeStep The number of the time step
 * @param iSurface The index of the surface position
 * @param nInter The quantity of interstitial at each surface position
 * @param previousFlux The previous I flux at each surface position
 */
void writeSurface1D(int timeStep, int iSurface, double nInter,
		double previousFlux);

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
void writeSurface3D(int timeStep, std::vector<std::vector<int> > iSurface,
		std::vector<std::vector<double> > nInter,
		std::vector<std::vector<double> > previousFlux);

/**
 * Write the bottom informations in the
 * concentration subgroup.
 *
 * @param nHe The quantity of helium at the bottom
 * @param previousHeFlux The previous He flux
 * @param nD The quantity of deuterium at the bottom
 * @param previousDFlux The previous D flux
 * @param nT The quantity of tritium at the bottom
 * @param previousTFlux The previous T flux
 */
void writeBottom1D(double nHe, double previousHeFlux, double nD,
		double previousDFlux, double nT, double previousTFlux);

/**
 * Write the bottom informations in the
 * concentration subgroup.
 *
 * @param nHe The quantity of helium at the bottom
 * @param previousHeFlux The previous He flux
 * @param nD The quantity of deuterium at the bottom
 * @param previousDFlux The previous D flux
 * @param nT The quantity of tritium at the bottom
 * @param previousTFlux The previous T flux
 */
void writeBottom2D(std::vector<double> nHe, std::vector<double> previousHeFlux,
		std::vector<double> nD, std::vector<double> previousDFlux,
		std::vector<double> nT, std::vector<double> previousTFlux);

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
void fillConcentrations(double concArray[][2], int i, int j = -1, int k = -1);

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
std::vector<std::vector<int> > readSurface3D(const std::string& fileName,
		int lastTimeStep);

/**
 * Read some data from the concentration
 * group of a HDF5 file in the case of a 1D grid (one float).
 *
 * @param fileName The name of the file to read from
 * @param lastTimeStep The value of the last written time step
 * @param dataName The name of the data we want
 * @return The value of the data
 */
double readData1D(const std::string& fileName, int lastTimeStep,
		const std::string& dataName);

/**
 * Read some data from the concentration group of a HDF5 file in
 * the case of a 2D grid (a vector).
 *
 * @param fileName The name of the file to read from
 * @param lastTimeStep The value of the last written time step
 * @param dataName The name of the data we want
 * @return The vector of the data
 */
std::vector<double> readData2D(const std::string& fileName, int lastTimeStep,
		const std::string& dataName);

/**
 * Read some data from the concentration group of a HDF5 file in
 * the case of a 3D grid (a vector of vector).
 *
 * @param fileName The name of the file to read from
 * @param lastTimeStep The value of the last written time step
 * @param dataName The name of the data we want
 * @return The vector of vector of data
 */
std::vector<std::vector<double> > readData3D(const std::string& fileName,
		int lastTimeStep, const std::string& dataName);

/**
 * Read the network sizes from a HDF5 file.
 * @param fileName The name of the file to read from.
 * @param normalSize The number of normal clusters.
 * @param superSize The number of super clusters.
 */
void readNetworkSize(const std::string& fileName, int &normalSize,
		int &superSize);

/**
 * Read the cluster from a HDF5 file.
 * @param id The id of the cluster.
 * @param formationEnergy The formation energy.
 * @param migrationEnergy The migration energy.
 * @param diffusionFactor The diffusion factor.
 * @return The cluster composition.
 */
std::vector<int> readCluster(int id, double &formationEnergy,
		double &migrationEnergy, double &diffusionFactor);

/**
 * Read the super cluster from a HDF5 file.
 * @param id The id of the cluster.
 * @return The list of clusters that it contains.
 */
std::set<std::tuple<int, int, int, int> > readSuperCluster(int id);

/**
 * Read and set the reactions in the network.
 *
 * @param network The network in which the reactions will be set
 */
void readReactions(IReactionNetwork& network);

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
std::vector<std::vector<double> > readGridPoint(const std::string& fileName,
		int lastTimeStep, int i, int j = -1, int k = -1);

}

} /* namespace xolotlCore */
#endif
