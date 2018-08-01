#include "HDF5Utils.h"
#include <PSICluster.h>
#include <PSISuperCluster.h>
#include <iostream>
#include <sstream>
#include <hdf5.h>
#include <mpi.h>

using namespace xolotlCore;

hid_t propertyListId, fileId, concentrationGroupId, subConcGroupId,
		concDataspaceId, headerGroupId, networkGroupId;
herr_t status;

void HDF5Utils::initializeFile(const std::string& fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Create the file
	fileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Create the group where the header will be stored
	headerGroupId = H5Gcreate2(fileId, "headerGroup", H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);

	// Create the group where the concentrations will be stored
	concentrationGroupId = H5Gcreate2(fileId, "concentrationsGroup",
	H5P_DEFAULT,
	H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the last written time step attribute
	int lastTimeStep = -1;
	hid_t lastDataspaceId = H5Screate(H5S_SCALAR);
	hid_t lastAttributeId = H5Acreate2(concentrationGroupId, "lastTimeStep",
	H5T_STD_I32LE, lastDataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);

	// Close everything
	status = H5Aclose(lastAttributeId);
	status = H5Sclose(lastDataspaceId);

	return;
}

void HDF5Utils::openFile(const std::string& fileName) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read and write access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the concentration group
	concentrationGroupId = H5Gopen(fileId, "concentrationsGroup", H5P_DEFAULT);
	return;
}

void HDF5Utils::fillHeader(std::vector<double>& grid, int ny, double hy, int nz,
		double hz) {
	// Create, write, and close the nx attribute
	int nx = grid.size() - 2;
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(headerGroupId, "nx", H5T_STD_I32LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Create, write, and close the hx attribute
	double hx = 0.0;
	if (grid.size() > 0)
		hx = grid[1] - grid[0];
	attributeId = H5Acreate2(headerGroupId, "hx", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Create, write, and close the ny attribute
	attributeId = H5Acreate2(headerGroupId, "ny", H5T_STD_I32LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Create, write, and close the hy attribute
	attributeId = H5Acreate2(headerGroupId, "hy", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Create, write, and close the nz attribute
	attributeId = H5Acreate2(headerGroupId, "nz", H5T_STD_I32LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Create, write, and close the hz attribute
	attributeId = H5Acreate2(headerGroupId, "hz", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hz);

	// Create, write, and close the grid dataset
	double gridArray[nx];
	for (int i = 0; i < nx; i++) {
		gridArray[i] = grid[i + 1] - grid[1];
	}
	hsize_t dims[1];
	dims[0] = nx;
	dataspaceId = H5Screate_simple(1, dims, NULL);
	hid_t datasetId = H5Dcreate2(headerGroupId, "grid", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &gridArray);
	status = H5Dclose(datasetId);

	// Close everything
	status = H5Aclose(attributeId);
	status = H5Sclose(dataspaceId);

	return;
}
void HDF5Utils::fillNetworkComp(std::vector<std::vector<int> > compVec) {
	// Create the array that will store the compositions and fill it
	int dof = compVec.size();
	int compSize = compVec[0].size();
	int compArray[dof][compSize];
	for (int i = 0; i < dof; i++) {
		for (int j = 0; j < compSize; j++) {
			compArray[i][j] = compVec[i][j];
		}
	}

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[2];
	dims[0] = dof;
	dims[1] = compSize;
	hid_t dataspaceId = H5Screate_simple(2, dims, NULL);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(headerGroupId, "composition", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &compArray);
	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::fillNetwork(const std::string& fileName,
		IReactionNetwork& network) {
	// Initial declaration
	bool groupExist = false;

	// Check if the group was already saved
	if (!fileName.empty()) {
		// Set up file access property list with parallel I/O access
		propertyListId = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

		// Open the given HDF5 file with read only access
		hid_t fromFileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY,
				propertyListId);

		// Close the property list
		status = H5Pclose(propertyListId);

		// Check the group
		groupExist = H5Lexists(fromFileId, "/networkGroup", H5P_DEFAULT);

		// If the group exist
		if (groupExist) {
			// Copy it
			status = H5Ocopy(fromFileId, "/networkGroup", fileId,
					"/networkGroup",
					H5P_DEFAULT, H5P_DEFAULT);
		}

		// Close the from file
		status = H5Fclose(fromFileId);
	} else if (!groupExist) {
		// Write it from scratch
		// Create the group and its size attributes
		networkGroupId = H5Gcreate2(fileId, "networkGroup", H5P_DEFAULT,
		H5P_DEFAULT,
		H5P_DEFAULT);
		int totalSize = network.size(), superSize = network.getAll(
				ReactantType::PSISuper).size(), normalSize = totalSize
				- superSize;
		hid_t dataspaceId = H5Screate(H5S_SCALAR);
		hid_t attributeId = H5Acreate2(networkGroupId, "normalSize",
		H5T_STD_I32LE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attributeId, H5T_STD_I32LE, &normalSize);
		status = H5Aclose(attributeId);
		attributeId = H5Acreate2(networkGroupId, "superSize", H5T_STD_I32LE,
				dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attributeId, H5T_STD_I32LE, &superSize);
		status = H5Aclose(attributeId);

		// Loop on all the clusters
		auto& allReactants = network.getAll();
		std::for_each(allReactants.begin(), allReactants.end(),
				[](IReactant& currReactant) {
					// Create a group named after its id
					int id = currReactant.getId() - 1;
					std::stringstream subGroupName;
					subGroupName << id;
					hid_t clusterGroupId = H5Gcreate2(networkGroupId, subGroupName.str().c_str(), H5P_DEFAULT, H5P_DEFAULT,
							H5P_DEFAULT);
					// Super cluster case
					if (currReactant.getType() == ReactantType::PSISuper) {
						// Write the dataset with the coordinate of each contained cluster
						auto& currCluster = static_cast<PSISuperCluster&>(currReactant);
						int nTot = currCluster.getNTot();
						int heVArray[nTot][4];
						auto& heVList = currCluster.getCoordList();
						int i = 0;
						for (auto const& pair : heVList) {
							heVArray[i][0] = std::get<0>(pair);
							heVArray[i][1] = std::get<1>(pair);
							heVArray[i][2] = std::get<2>(pair);
							heVArray[i][3] = std::get<3>(pair);
							i++;
						}
						hsize_t dim[2]; dim[0] = nTot; dim[1] = 4;
						hid_t dataspaceId = H5Screate_simple(2, dim, NULL);
						hid_t datasetId = H5Dcreate2(clusterGroupId, "heVList", H5T_STD_I32LE,
								dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
								H5P_DEFAULT, &heVArray);
						status = H5Dclose(datasetId);
						status = H5Sclose(dataspaceId);
					}
					// Normal cluster case
					else {
						// Write the composition attribute
						auto& comp = currReactant.getComposition();
						hsize_t dim[1]; dim[0] = comp.size();
						hid_t dataspaceId = H5Screate_simple(1, dim, NULL);
						int compArray[dim[0]];
						for (int i = 0; i < dim[0]; i++) {
							compArray[i] = comp[i];
						}
						hid_t attributeId = H5Acreate2(clusterGroupId, "composition",
								H5T_STD_I32LE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
						status = H5Awrite(attributeId, H5T_STD_I32LE, &compArray);
						status = H5Aclose(attributeId);
						status = H5Sclose(dataspaceId);

						// Write the energy and diffusion parameters
						dataspaceId = H5Screate(H5S_SCALAR);
						attributeId = H5Acreate2(clusterGroupId, "formationEnergy",
								H5T_IEEE_F64LE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
						double formationEnergy = currReactant.getFormationEnergy();
						status = H5Awrite(attributeId, H5T_IEEE_F64LE, &formationEnergy);
						status = H5Aclose(attributeId);
						attributeId = H5Acreate2(clusterGroupId, "migrationEnergy",
								H5T_IEEE_F64LE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
						double migrationEnergy = currReactant.getMigrationEnergy();
						status = H5Awrite(attributeId, H5T_IEEE_F64LE, &migrationEnergy);
						status = H5Aclose(attributeId);
						attributeId = H5Acreate2(clusterGroupId, "diffusionFactor",
								H5T_IEEE_F64LE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
						double diffusionFactor = currReactant.getDiffusionFactor();
						status = H5Awrite(attributeId, H5T_IEEE_F64LE, &diffusionFactor);
						status = H5Aclose(attributeId);
						status = H5Sclose(dataspaceId);
					}

					// Write a dataset of production reactions
					auto prodVec = currReactant.getProdVector();
					if (prodVec.size() > 0) {
						int nReact = prodVec.size();
						int dataSize = prodVec[0].size();
						double prodArray[nReact][dataSize];
						for (int i = 0; i < nReact; i++) {
							for (int j = 0; j < dataSize; j++) {
								prodArray[i][j] = prodVec[i][j];
							}
						}
						hsize_t dims[2];
						dims[0] = nReact;
						dims[1] = dataSize;
						hid_t dataspaceId = H5Screate_simple(2, dims, NULL);
						hid_t datasetId = H5Dcreate2(clusterGroupId, "prod", H5T_IEEE_F64LE,
								dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
								H5P_DEFAULT, &prodArray);
						status = H5Dclose(datasetId);
						status = H5Sclose(dataspaceId);
					}

					// Write a dataset of combination reactions
					auto combVec = currReactant.getCombVector();
					if (combVec.size() > 0) {
						int nReact = combVec.size();
						int dataSize = combVec[0].size();
						double combArray[nReact][dataSize];
						for (int i = 0; i < nReact; i++) {
							for (int j = 0; j < dataSize; j++) {
								combArray[i][j] = combVec[i][j];
							}
						}
						hsize_t dims[2];
						dims[0] = nReact;
						dims[1] = dataSize;
						hid_t dataspaceId = H5Screate_simple(2, dims, NULL);
						hid_t datasetId = H5Dcreate2(clusterGroupId, "comb", H5T_IEEE_F64LE,
								dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
								H5P_DEFAULT, &combArray);
						status = H5Dclose(datasetId);
						status = H5Sclose(dataspaceId);
					}

					// Write a dataset of dissociation reactions
					auto dissoVec = currReactant.getDissoVector();
					if (dissoVec.size() > 0) {
						int nReact = dissoVec.size();
						int dataSize = dissoVec[0].size();
						double dissoArray[nReact][dataSize];
						for (int i = 0; i < nReact; i++) {
							for (int j = 0; j < dataSize; j++) {
								dissoArray[i][j] = dissoVec[i][j];
							}
						}
						hsize_t dims[2];
						dims[0] = nReact;
						dims[1] = dataSize;
						hid_t dataspaceId = H5Screate_simple(2, dims, NULL);
						hid_t datasetId = H5Dcreate2(clusterGroupId, "disso", H5T_IEEE_F64LE,
								dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
								H5P_DEFAULT, &dissoArray);
						status = H5Dclose(datasetId);
						status = H5Sclose(dataspaceId);
					}

					// Write a dataset of emission reactions
					auto emitVec = currReactant.getEmitVector();
					if (emitVec.size() > 0) {
						int nReact = emitVec.size();
						int dataSize = emitVec[0].size();
						double emitArray[nReact][dataSize];
						for (int i = 0; i < nReact; i++) {
							for (int j = 0; j < dataSize; j++) {
								emitArray[i][j] = emitVec[i][j];
							}
						}
						hsize_t dims[2];
						dims[0] = nReact;
						dims[1] = dataSize;
						hid_t dataspaceId = H5Screate_simple(2, dims, NULL);
						hid_t datasetId = H5Dcreate2(clusterGroupId, "emit", H5T_IEEE_F64LE,
								dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
								H5P_DEFAULT, &emitArray);
						status = H5Dclose(datasetId);
						status = H5Sclose(dataspaceId);
					}

					// Close the group
					status = H5Gclose(clusterGroupId);
				});

		// Close everything
		status = H5Sclose(dataspaceId);
		status = H5Gclose(networkGroupId);
	}

	return;
}

void HDF5Utils::addConcentrationSubGroup(int timeStep, double time,
		double previousTime, double deltaTime) {
	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentration_" << timeStep;

	// Create the subgroup where the concentrations at this time step will be stored
	subConcGroupId = H5Gcreate2(concentrationGroupId,
			subGroupName.str().c_str(),
			H5P_DEFAULT,
			H5P_DEFAULT, H5P_DEFAULT);

	// Create, write, and close the absolute time attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(subConcGroupId, "absoluteTime",
	H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(attributeId);

	// Create, write, and close the previous time attribute
	attributeId = H5Acreate2(subConcGroupId, "previousTime", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousTime);
	status = H5Aclose(attributeId);

	// Create, write, and close the timestep time attribute
	attributeId = H5Acreate2(subConcGroupId, "deltaTime", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(attributeId);
	status = H5Sclose(dataspaceId);

	// Overwrite the last time step attribute of the concentration group
	attributeId = H5Aopen(concentrationGroupId, "lastTimeStep", H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &timeStep);
	status = H5Aclose(attributeId);

	return;
}

void HDF5Utils::writeSurface1D(int timeStep, int iSurface, double nInter,
		double previousFlux) {
	// Create, write, and close the surface position attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(subConcGroupId, "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &iSurface);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of interstitial attribute
	attributeId = H5Acreate2(subConcGroupId, "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nInter);
	status = H5Aclose(attributeId);

	// Create, write, and close the flux of interstitial attribute
	attributeId = H5Acreate2(subConcGroupId, "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousFlux);
	status = H5Aclose(attributeId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::writeSurface2D(int timeStep, std::vector<int> iSurface,
		std::vector<double> nInter, std::vector<double> previousFlux) {
	// Create the array that will store the indices and fill it
	int size = iSurface.size();
	int indexArray[size];
	for (int i = 0; i < size; i++) {
		indexArray[i] = iSurface[i];
	}

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[1];
	dims[0] = size;
	hid_t dataspaceId = H5Screate_simple(1, dims, NULL);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(subConcGroupId, "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write surface array in the dataset
	status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &indexArray);

	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the array that will store the quantities and fill it
	double quantityArray[size];
	for (int i = 0; i < size; i++) {
		quantityArray[i] = nInter[i];
	}

	// Create the dataset for the surface indices
	datasetId = H5Dcreate2(subConcGroupId, "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write quantityArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);

	// Close the dataset
	status = H5Dclose(datasetId);

	// Fill the array with the previous flux
	for (int i = 0; i < size; i++) {
		quantityArray[i] = previousFlux[i];
	}

	// Create the dataset for the surface indices
	datasetId = H5Dcreate2(subConcGroupId, "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write quantityArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::writeSurface3D(int timeStep,
		std::vector<std::vector<int> > iSurface,
		std::vector<std::vector<double> > nInter,
		std::vector<std::vector<double> > previousFlux) {
	// Create the array that will store the indices and fill it
	int xSize = iSurface.size();
	int ySize = iSurface[0].size();
	int indexArray[xSize][ySize];
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			indexArray[i][j] = iSurface[i][j];
		}
	}

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[2];
	dims[0] = xSize;
	dims[1] = ySize;
	hid_t dataspaceId = H5Screate_simple(2, dims, NULL);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(subConcGroupId, "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &indexArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the array that will store the interstitial quantities and fill it
	double quantityArray[xSize][ySize];
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			quantityArray[i][j] = nInter[i][j];
		}
	}

	// Create the dataset for the interstitial quantities
	datasetId = H5Dcreate2(subConcGroupId, "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Fill the array that will store the interstitial flux
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			quantityArray[i][j] = previousFlux[i][j];
		}
	}

	// Create the dataset for the interstitial quantities
	datasetId = H5Dcreate2(subConcGroupId, "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::writeBottom1D(double nHe, double previousHeFlux, double nD,
		double previousDFlux, double nT, double previousTFlux) {
	// Create, write, and close the quantity of helium attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(subConcGroupId, "nHelium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nHe);
	status = H5Aclose(attributeId);
	// Create, write, and close the flux of helium attribute
	attributeId = H5Acreate2(subConcGroupId, "previousHeFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousHeFlux);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of deuterium attribute
	attributeId = H5Acreate2(subConcGroupId, "nDeuterium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nD);
	status = H5Aclose(attributeId);
	// Create, write, and close the flux of deuterium attribute
	attributeId = H5Acreate2(subConcGroupId, "previousDFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousDFlux);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of tritium attribute
	attributeId = H5Acreate2(subConcGroupId, "nTritium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nT);
	status = H5Aclose(attributeId);
	// Create, write, and close the flux of tritium attribute
	attributeId = H5Acreate2(subConcGroupId, "previousTFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousTFlux);
	status = H5Aclose(attributeId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::writeBottom2D(std::vector<double> nHe,
		std::vector<double> previousHeFlux, std::vector<double> nD,
		std::vector<double> previousDFlux, std::vector<double> nT,
		std::vector<double> previousTFlux) {
	// Find out the size of the arrays
	const int size = nHe.size();

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[1];
	dims[0] = size;
	hid_t dataspaceId = H5Screate_simple(1, dims, NULL);

	// Create the dataset for helium
	hid_t datasetId = H5Dcreate2(subConcGroupId, "nHelium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	// Create the array that will store the quantities and fill it
//	double quantityArray[size];
//	for (int i = 0; i < size; i++) {
//		quantityArray[i] = nInter[i];
//	}
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nHe.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the helium flux
	datasetId = H5Dcreate2(subConcGroupId, "previousHeFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousHeFlux.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the deuterium
	datasetId = H5Dcreate2(subConcGroupId, "nDeuterium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nD.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the deuterium flux
	datasetId = H5Dcreate2(subConcGroupId, "previousDFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousDFlux.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the tritium
	datasetId = H5Dcreate2(subConcGroupId, "nTritium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nT.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the tritium flux
	datasetId = H5Dcreate2(subConcGroupId, "previousTFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousTFlux.data());

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);

	return;
}

void HDF5Utils::addConcentrationDataset(int size, int i, int j, int k) {
	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Create the dataspace for the dataset with dimension dims
	hsize_t dims[2];
	dims[0] = size;
	dims[1] = 2;
	concDataspaceId = H5Screate_simple(2, dims, NULL);

	// Create the dataset of concentrations for this position
	hid_t datasetId = H5Dcreate2(subConcGroupId, datasetName.str().c_str(),
	H5T_IEEE_F64LE, concDataspaceId,
	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Close everything
	status = H5Sclose(concDataspaceId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::fillConcentrations(double concArray[][2], int i, int j, int k) {
	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Open the already created dataset of concentrations for this position
	hid_t datasetId = H5Dopen(subConcGroupId, datasetName.str().c_str(),
	H5P_DEFAULT);

	// Write concArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, concArray);

	// Close everything
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::finalizeFile() {
	// Close everything
	status = H5Gclose(headerGroupId);
	status = H5Gclose(concentrationGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::closeFile() {
	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Gclose(concentrationGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::readHeader(const std::string& fileName, int &nx, double &hx,
		int &ny, double &hy, int &nz, double &hz) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the header group
	hid_t groupId = H5Gopen(fileId, "/headerGroup", H5P_DEFAULT);

	// Open and read the nx attribute
	hid_t attributeId = H5Aopen(groupId, "nx", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Open and read the hx attribute
	attributeId = H5Aopen(groupId, "hx", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Open and read the ny attribute
	attributeId = H5Aopen(groupId, "ny", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Open and read the hy attribute
	attributeId = H5Aopen(groupId, "hy", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Open and read the nz attribute
	attributeId = H5Aopen(groupId, "nz", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Open and read the hz attribute
	attributeId = H5Aopen(groupId, "hz", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hz);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(groupId);
	status = H5Fclose(fileId);

	return;
}

bool HDF5Utils::hasConcentrationGroup(const std::string& fileName,
		int &lastTimeStep) {
	// Initialize the boolean to return
	bool hasGroup = true;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Check the group
	bool groupExist = H5Lexists(fileId, "/concentrationsGroup", H5P_DEFAULT);
	// If the group exist
	if (groupExist) {
		// Open the concentration group
		concentrationGroupId = H5Gopen(fileId, "/concentrationsGroup",
		H5P_DEFAULT);

		// Open and read the lastTimeStep attribute
		hid_t lastAttributeId = H5Aopen(concentrationGroupId, "lastTimeStep",
		H5P_DEFAULT);
		status = H5Aread(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);
		status = H5Aclose(lastAttributeId);

		status = H5Gclose(concentrationGroupId);

		// if lastTimeStep is still negative the group is not valid
		if (lastTimeStep < 0)
			hasGroup = false;
	}
	// if not
	else
		hasGroup = false;

	// Close everything
	status = H5Fclose(fileId);

	return hasGroup;
}

void HDF5Utils::readTimes(const std::string& fileName, int lastTimeStep,
		double &time, double &deltaTime) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open and read the absoluteTime attribute
	hid_t attributeId = H5Aopen(subConcGroupId, "absoluteTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(attributeId);

	// Open and read the deltaTime attribute
	attributeId = H5Aopen(subConcGroupId, "deltaTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return;
}

double HDF5Utils::readPreviousTime(const std::string& fileName,
		int lastTimeStep) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open and read the previousTime attribute
	double previousTime = 0.0;
	hid_t attributeId = H5Aopen(subConcGroupId, "previousTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &previousTime);

	// Close everything
	status = H5Aclose(attributeId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return previousTime;
}

int HDF5Utils::readSurface1D(const std::string& fileName, int lastTimeStep) {
	// Initialize the surface position
	int iSurface = 0;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(subConcGroupId, "iSurface", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &iSurface);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return iSurface;
}

std::vector<int> HDF5Utils::readSurface2D(const std::string& fileName,
		int lastTimeStep) {
	// Create the vector to return
	std::vector<int> toReturn;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "iSurface", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[1];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	int index[dims[0]];

	// Read the data set
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &index);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		toReturn.push_back(index[i]);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

std::vector<std::vector<int> > HDF5Utils::readSurface3D(
		const std::string& fileName, int lastTimeStep) {
	// Create the vector to return
	std::vector<std::vector<int> > toReturn;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, "iSurface", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[2];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	int index[dims[0]][dims[1]];

	// Read the data set
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &index);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		// Create a temporary vector
		std::vector<int> temp;
		for (int j = 0; j < dims[1]; j++) {
			temp.push_back(index[i][j]);
		}
		// Add the temporary vector to the one to return
		toReturn.push_back(temp);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

double HDF5Utils::readData1D(const std::string& fileName, int lastTimeStep,
		const std::string& dataName) {
	// Initialize the surface position
	double data = 0.0;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(subConcGroupId, dataName.c_str(), H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &data);
	status = H5Aclose(attributeId);

	// Close everything
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return data;
}

std::vector<double> HDF5Utils::readData2D(const std::string& fileName,
		int lastTimeStep, const std::string& dataName) {
	// Create the vector to return
	std::vector<double> toReturn;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, dataName.c_str(), H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[1];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	double flux[dims[0]];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &flux);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		toReturn.push_back(flux[i]);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

std::vector<std::vector<double> > HDF5Utils::readData3D(
		const std::string& fileName, int lastTimeStep,
		const std::string& dataName) {
	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the name of the sub group
	std::stringstream subGroupName;
	subGroupName << "concentrationsGroup/concentration_" << lastTimeStep;

	// Open this specific concentration sub group
	subConcGroupId = H5Gopen(fileId, subGroupName.str().c_str(), H5P_DEFAULT);

	// Open the dataset
	hid_t datasetId = H5Dopen(subConcGroupId, dataName.c_str(), H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	hsize_t dims[2];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

	// Create the array that will receive the indices
	double quantity[dims[0]][dims[1]];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantity);

	// Loop on the length and fill the vector to return
	for (int i = 0; i < dims[0]; i++) {
		// Create a temporary vector
		std::vector<double> temp;
		for (int j = 0; j < dims[1]; j++) {
			temp.push_back(quantity[i][j]);
		}
		// Add the temporary vector to the one to return
		toReturn.push_back(temp);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
	status = H5Gclose(subConcGroupId);
	status = H5Fclose(fileId);

	return toReturn;
}

void HDF5Utils::readNetworkSize(const std::string& fileName, int &normalSize,
		int &superSize) {
	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Open the network group
	networkGroupId = H5Gopen(fileId, "networkGroup", H5P_DEFAULT);

	// Open and read the normalSize attribute
	hid_t attributeId = H5Aopen(networkGroupId, "normalSize", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &normalSize);
	status = H5Aclose(attributeId);

	// Open and read the superSize attribute
	attributeId = H5Aopen(networkGroupId, "superSize", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &superSize);
	status = H5Aclose(attributeId);

	return;
}

std::vector<int> HDF5Utils::readCluster(int id, double &formationEnergy,
		double &migrationEnergy, double &diffusionFactor) {
	// Open the corresponding group
	std::stringstream idName;
	idName << id;
	hid_t clusterGroupId = H5Gopen(networkGroupId, idName.str().c_str(),
	H5P_DEFAULT);

	// Open and read the formation energy attribute
	hid_t attributeId = H5Aopen(clusterGroupId, "formationEnergy", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &formationEnergy);
	status = H5Aclose(attributeId);

	// Open and read the migration energy attribute
	attributeId = H5Aopen(clusterGroupId, "migrationEnergy", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &migrationEnergy);
	status = H5Aclose(attributeId);

	// Open and read the diffusion factor attribute
	attributeId = H5Aopen(clusterGroupId, "diffusionFactor", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &diffusionFactor);
	status = H5Aclose(attributeId);

	// Read the composition attribute
	attributeId = H5Aopen_name(clusterGroupId, "composition");
	hid_t dataspaceId = H5Aget_space(attributeId);
	hsize_t dims[1];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
	int compArray[dims[0]];
	status = H5Aread(attributeId, H5T_STD_I32LE, &compArray);
	status = H5Aclose(attributeId);

	// Fill the vector
	std::vector<int> comp;
	for (int i = 0; i < dims[0]; i++) {
		comp.push_back(compArray[i]);
	}

	return comp;
}

std::set<std::tuple<int, int, int, int> > HDF5Utils::readSuperCluster(int id) {
	// Open the corresponding group
	std::stringstream idName;
	idName << id;
	hid_t clusterGroupId = H5Gopen(networkGroupId, idName.str().c_str(),
	H5P_DEFAULT);

	// Read the heVList dataset
	hid_t datasetId = H5Dopen(clusterGroupId, "heVList", H5P_DEFAULT);
	hid_t dataspaceId = H5Dget_space(datasetId);
	hsize_t dims[2];
	status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
	int heVArray[dims[0]][dims[1]];
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &heVArray);
	status = H5Dclose(datasetId);

	// Fill the set
	std::set<std::tuple<int, int, int, int> > heVList;
	for (int i = 0; i < dims[0]; i++) {
		auto pair = std::make_tuple(heVArray[i][0], heVArray[i][1],
				heVArray[i][2], heVArray[i][3]);
		heVList.emplace(pair);
	}

	return heVList;
}

void HDF5Utils::readReactions(IReactionNetwork& network) {
	// Loop on the reactants
	auto& allReactants = network.getAll();
	std::for_each(allReactants.begin(), allReactants.end(),
			[&allReactants, &network](IReactant& currReactant) {
				// Open the corresponding group
				int id = currReactant.getId() - 1;
				std::stringstream idName;
				idName << id;
				hid_t clusterGroupId = H5Gopen(networkGroupId, idName.str().c_str(),
						H5P_DEFAULT);

				// Read the production dataset
				bool datasetExist = H5Lexists(clusterGroupId, "prod",
						H5P_DEFAULT);
				if (datasetExist) {
					hid_t datasetId = H5Dopen(clusterGroupId, "prod", H5P_DEFAULT);
					hid_t dataspaceId = H5Dget_space(datasetId);
					hsize_t dims[2];
					status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
					double prodVec[dims[0]][dims[1]];
					status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
							H5P_DEFAULT, &prodVec);
					status = H5Dclose(datasetId);
					// Loop on the prod vector
					for (int i = 0; i < dims[0]; i++) {
						// Get pointers to the 2 reactants
						auto& firstReactant = allReactants.at(prodVec[i][0]);
						auto& secondReactant = allReactants.at(prodVec[i][1]);

						// Create and add the reaction to the network
						std::unique_ptr<ProductionReaction> reaction(
								new ProductionReaction(firstReactant, secondReactant));
						auto& prref = network.add(std::move(reaction));

						// Add the reaction to the cluster
						currReactant.resultFrom(prref, prodVec[i][2]);
					}
				}

				// Read the combination dataset
				datasetExist = H5Lexists(clusterGroupId, "comb",
						H5P_DEFAULT);
				if (datasetExist) {
					hid_t datasetId = H5Dopen(clusterGroupId, "comb", H5P_DEFAULT);
					hid_t dataspaceId = H5Dget_space(datasetId);
					hsize_t dims[2];
					status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
					double combVec[dims[0]][dims[1]];
					status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
							H5P_DEFAULT, &combVec);
					status = H5Dclose(datasetId);
					// Loop on the prod vector
					for (int i = 0; i < dims[0]; i++) {
						// Get pointers to the combining reactant
						auto& firstReactant = allReactants.at(combVec[i][0]);

						// Create and add the reaction to the network
						std::unique_ptr<ProductionReaction> reaction(
								new ProductionReaction(firstReactant, currReactant));
						auto& prref = network.add(std::move(reaction));

						// Add the reaction to the cluster
						currReactant.participateIn(prref, combVec[i][1]);
					}
				}

				// Read the dissociation dataset
				datasetExist = H5Lexists(clusterGroupId, "disso",
						H5P_DEFAULT);
				if (datasetExist) {
					hid_t datasetId = H5Dopen(clusterGroupId, "disso", H5P_DEFAULT);
					hid_t dataspaceId = H5Dget_space(datasetId);
					hsize_t dims[2];
					status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
					double dissoVec[dims[0]][dims[1]];
					status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
							H5P_DEFAULT, &dissoVec);
					status = H5Dclose(datasetId);
					// Loop on the prod vector
					for (int i = 0; i < dims[0]; i++) {
						// Get pointers to the other reactants
						auto& emittingReactant = allReactants.at(dissoVec[i][0]);
						auto& secondReactant = allReactants.at(dissoVec[i][1]);

						// Create and add the reaction to the network
						std::unique_ptr<ProductionReaction> reaction(
								new ProductionReaction(currReactant, secondReactant));
						auto& prref = network.add(std::move(reaction));
						std::unique_ptr<DissociationReaction> dissociationReaction(
								new DissociationReaction(emittingReactant, prref.first,
										prref.second, &prref));
						auto& drref = network.add(std::move(dissociationReaction));

						// Add the reaction to the cluster
						currReactant.participateIn(drref, dissoVec[i][2]);
					}
				}

				// Read the emission dataset
				datasetExist = H5Lexists(clusterGroupId, "emit",
						H5P_DEFAULT);
				if (datasetExist) {
					hid_t datasetId = H5Dopen(clusterGroupId, "emit", H5P_DEFAULT);
					hid_t dataspaceId = H5Dget_space(datasetId);
					hsize_t dims[2];
					status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
					double emitVec[dims[0]][dims[1]];
					status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
							H5P_DEFAULT, &emitVec);
					status = H5Dclose(datasetId);
					// Loop on the prod vector
					for (int i = 0; i < dims[0]; i++) {
						// Get pointers to the other reactants
						auto& firstReactant = allReactants.at(emitVec[i][0]);
						auto& secondReactant = allReactants.at(emitVec[i][1]);

						// Create and add the reaction to the network
						std::unique_ptr<ProductionReaction> reaction(
								new ProductionReaction(firstReactant, secondReactant));
						auto& prref = network.add(std::move(reaction));
						std::unique_ptr<DissociationReaction> dissociationReaction(
								new DissociationReaction(currReactant, prref.first,
										prref.second, &prref));
						auto& drref = network.add(std::move(dissociationReaction));

						// Add the reaction to the cluster
						currReactant.emitFrom(drref, emitVec[i][2]);
					}
				}
			});
}

std::vector<std::vector<double> > HDF5Utils::readGridPoint(
		const std::string& fileName, int lastTimeStep, int i, int j, int k) {
	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

	// Set up file access property list with parallel I/O access
	propertyListId = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(propertyListId, MPI_COMM_WORLD, MPI_INFO_NULL);

	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, propertyListId);

	// Close the property list
	status = H5Pclose(propertyListId);

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "concentrationsGroup/concentration_" << lastTimeStep
			<< "/position_" << i << "_" << j << "_" << k;

	// Check the dataset
	bool datasetExist = H5Lexists(fileId, datasetName.str().c_str(),
	H5P_DEFAULT);

	// If the dataset exists
	if (datasetExist) {
		// Open the dataset
		hid_t datasetId = H5Dopen(fileId, datasetName.str().c_str(),
		H5P_DEFAULT);

		// Get the dataspace object
		hid_t dataspaceId = H5Dget_space(datasetId);

		// Get the dimensions of the dataset
		hsize_t dims[2];
		status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);

		// Create the array that will receive the concentrations
		double conc[dims[0]][dims[1]];

		// Read the data set
		status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, &conc);

		// Loop on the length
		for (unsigned int n = 0; n < dims[0]; n++) {
			// Create the concentration vector for this cluster
			std::vector<double> tmp;
			tmp.push_back(conc[n][0]);
			tmp.push_back(conc[n][1]);

			// Add it to the main vector
			toReturn.push_back(tmp);
		}

		// Close everything
		status = H5Dclose(datasetId);
		status = H5Sclose(dataspaceId);
	}

	// Close the file
	status = H5Fclose(fileId);

	return toReturn;
}
