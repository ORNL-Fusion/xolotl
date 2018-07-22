#include <iostream>
#include <sstream>
#include <iterator>
#include <array>
#include "hdf5.h"
#include "mpi.h"
#include "xolotlCore/io/XFile.h"

namespace xolotlCore {


HDF5File::AccessMode XFile::EnsureCreateAccessMode(HDF5File::AccessMode mode) {

    bool createMode = 
        ((mode == HDF5File::AccessMode::CreateOrTruncateIfExists) or 
            (mode == HDF5File::AccessMode::CreateOrFailIfExists));

    if(not createMode) {
        throw HDF5Exception("Attempt to create existing file with a non-create access mode");
    }
    return mode;
}

HDF5File::AccessMode XFile::EnsureOpenAccessMode(HDF5File::AccessMode mode) {

    bool openMode = 
        ((mode == HDF5File::AccessMode::OpenReadOnly) or 
            (mode == HDF5File::AccessMode::OpenReadWrite));

    if(not openMode) {
        throw HDF5Exception("Attempt to open existing file with a non-open access mode");
    }
    return mode;
}


XFile::XFile(fs::path _path,
        const std::vector<double>& grid,
        const XFile::HeaderGroup::NetworkCompsType& compVec,
        fs::path networkFilePath,
        int ny,
        double hy,
        int nz,
        double hz,
        AccessMode _mode)
  : HDF5File(_path, EnsureCreateAccessMode(_mode), true)
{
    // Create and initialize the header group.
    HeaderGroup headerGroup(*this, grid, ny, hy, nz, hz, compVec);

	// Create and initialize the group where the concentrations will be stored
    ConcentrationGroup concGroup(*this, true);

    // Copy the network from another file if desired.
    if(not networkFilePath.empty()) {

        // Open the given file.
        XFile networkFile(networkFilePath, AccessMode::OpenReadOnly);

        // Check if given file has a network group.
        auto srcNetworkGroup = networkFile.getGroup<NetworkGroup>();
        if(srcNetworkGroup) {
            // Given file has a network group.  Copy it.
            srcNetworkGroup->copyTo(*this);
        }
    }
}


XFile::XFile(fs::path _path, AccessMode _mode)
  : HDF5File(_path, EnsureOpenAccessMode(_mode), true) {

    // Nothing else to do.
}


//----------------------------------------------------------------------------
// HeaderGroup
//
const fs::path XFile::HeaderGroup::path = "/headerGroup";
const std::string XFile::HeaderGroup::netCompsDatasetName = "composition";


XFile::HeaderGroup::HeaderGroup(const XFile& file,
        const std::vector<double>& grid,
        int ny, double hy, 
        int nz, double hz,
        const NetworkCompsType& compVec) 
  : HDF5File::Group(file, HeaderGroup::path, true) {

    // Base class created the group.

	// Create, write, and close the nx attribute
	int nx = grid.size() - 2;
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(getId(), "nx", H5T_STD_I32LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	auto status = H5Awrite(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Create, write, and close the hx attribute
	double hx = 0.0;
	if (grid.size() > 0)
		hx = grid[1] - grid[0];
	attributeId = H5Acreate2(getId(), "hx", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Create, write, and close the ny attribute
	attributeId = H5Acreate2(getId(), "ny", H5T_STD_I32LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Create, write, and close the hy attribute
	attributeId = H5Acreate2(getId(), "hy", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Create, write, and close the nz attribute
	attributeId = H5Acreate2(getId(), "nz", H5T_STD_I32LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Create, write, and close the hz attribute
	attributeId = H5Acreate2(getId(), "hz", H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &hz);

	// Create, write, and close the grid dataset
	double gridArray[nx];
	for (int i = 0; i < nx; i++) {
		gridArray[i] = grid[i + 1] - grid[1];
	}
    std::array<hsize_t, 1> dims{ (hsize_t)nx };
	dataspaceId = H5Screate_simple(dims.size(), dims.data(), nullptr);
	hid_t datasetId = H5Dcreate2(getId(), "grid", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &gridArray);
	status = H5Dclose(datasetId);

    // Initialize the network composition list.  Done here because
    // it is a dataset in the header group.
    initNetworkComps(compVec);

	// Close everything
	status = H5Aclose(attributeId);
	status = H5Sclose(dataspaceId);
}

XFile::HeaderGroup::HeaderGroup(const XFile& file)
  : HDF5File::Group(file, HeaderGroup::path, false) {

    // Base class opened the group, so nothing else to do.
}

void XFile::HeaderGroup::initNetworkComps(const NetworkCompsType& compVec) const {

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
    std::array<hsize_t, 2> dims{ (hsize_t)dof, (hsize_t)compSize };
	hid_t dataspaceId = H5Screate_simple(dims.size(), dims.data(), nullptr);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(getId(), netCompsDatasetName.c_str(), H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	auto  status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &compArray);
	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
}

void XFile::HeaderGroup::read(int &nx, double &hx,
		int &ny, double &hy, int &nz, double &hz) const {

	// Open and read the nx attribute
	hid_t attributeId = H5Aopen(getId(), "nx", H5P_DEFAULT);
	auto status = H5Aread(attributeId, H5T_STD_I32LE, &nx);
	status = H5Aclose(attributeId);
	// Open and read the hx attribute
	attributeId = H5Aopen(getId(), "hx", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hx);
	status = H5Aclose(attributeId);

	// Open and read the ny attribute
	attributeId = H5Aopen(getId(), "ny", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &ny);
	status = H5Aclose(attributeId);
	// Open and read the hy attribute
	attributeId = H5Aopen(getId(), "hy", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hy);
	status = H5Aclose(attributeId);

	// Open and read the nz attribute
	attributeId = H5Aopen(getId(), "nz", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_STD_I32LE, &nz);
	status = H5Aclose(attributeId);
	// Open and read the hz attribute
	attributeId = H5Aopen(getId(), "hz", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &hz);
	status = H5Aclose(attributeId);
}

XFile::HeaderGroup::NetworkCompsType
XFile::HeaderGroup::readNetworkComps(void) const {

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), netCompsDatasetName.c_str(), H5P_DEFAULT);

	// Get the dimensions of the dataset
    std::array<hsize_t, 2> dims;
	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);


	// Read the data set
    NetworkCompsType::value_type networkComps1D(dims[0] * dims[1]);
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			networkComps1D.data());

    // Convert the 1D dataset into 2D representation.
    NetworkCompsType networkComps(dims[0]);
    for(auto i = 0; i < dims[0]; ++i) {
        networkComps[i].reserve(dims[1]);
        std::copy(networkComps1D.begin() + i*dims[1],
                    networkComps1D.begin() + (i+1)*dims[1],
                    std::back_inserter(networkComps[i]));
    }

	// Close everything
	status = H5Dclose(datasetId);

	return networkComps;
}


//----------------------------------------------------------------------------
// NetworkGroup
//
const fs::path XFile::NetworkGroup::path = "/networkGroup";
const std::string XFile::NetworkGroup::netDatasetName = "network";


XFile::NetworkGroup::NetworkGroup(const XFile& file)
  : HDF5File::Group(file, NetworkGroup::path, false) {

    // Base class opened the group, so nothing else to do.
}


XFile::NetworkGroup::NetworkType XFile::NetworkGroup::readNetwork(void) const {

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), netDatasetName.c_str(), H5P_DEFAULT);

	// Get the dimensions of the dataset
    std::array<hsize_t, 2> dims;

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

	// Read the data set
    NetworkType::value_type network1D(dims[0] * dims[1]);
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			network1D.data());

	// Fill the vector to return with the dataset
    NetworkType network(dims[0]);
	// Loop on the size of the network
	for (auto i = 0; i < dims[0]; ++i) {
        network[i].reserve(dims[1]);
        std::copy(network1D.begin() + i*dims[1],
                    network1D.begin() + (i+1)*dims[1],
                    std::back_inserter(network[i]));
	}

	// Close everything
	status = H5Dclose(datasetId);

	return network;
}

void XFile::NetworkGroup::copyTo(const XFile& target) const {

    H5Ocopy(getLocation().getId(), NetworkGroup::path.string().c_str(),
            target.getId(), NetworkGroup::path.string().c_str(),
            H5P_DEFAULT,
            H5P_DEFAULT);
}



//----------------------------------------------------------------------------
// ConcentrationGroup
//
const fs::path XFile::ConcentrationGroup::path = "/concentrationsGroup";


XFile::ConcentrationGroup::ConcentrationGroup(const XFile& file, bool create)
  : HDF5File::Group(file, ConcentrationGroup::path, create) {

    if(create) {
        // We just created the group.
        // Populate it.
        //
        // Create, write, and close the last written time step attribute
        int lastTimeStep = -1;
        hid_t lastDataspaceId = H5Screate(H5S_SCALAR);
        hid_t lastAttributeId = H5Acreate2(getId(), "lastTimeStep",
        H5T_STD_I32LE, lastDataspaceId,
        H5P_DEFAULT, H5P_DEFAULT);
        auto status = H5Awrite(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);

        // Close everything
        status = H5Aclose(lastAttributeId);
        status = H5Sclose(lastDataspaceId);
    }
}


std::unique_ptr<XFile::TimestepGroup> 
XFile::ConcentrationGroup::addTimestepGroup(int timeStep,
                                        double time,
		                                double previousTime,
                                        double deltaTime) const {

    // Create a group for the new timestep.
    std::unique_ptr<XFile::TimestepGroup> tsGroup(
            new TimestepGroup(*this, timeStep, time, previousTime, deltaTime));

    // Update our last known timestep.
	auto attributeId = H5Aopen(getId(), "lastTimeStep", H5P_DEFAULT);
	auto status = H5Awrite(attributeId, H5T_STD_I32LE, &timeStep);
	status = H5Aclose(attributeId);

    return std::move(tsGroup);
}


int XFile::ConcentrationGroup::getLastTimeStep(void) const {

    int lastTimeStep;

    hid_t lastAttributeId = H5Aopen(getId(), "lastTimeStep", H5P_DEFAULT);
    auto status = H5Aread(lastAttributeId, H5T_STD_I32LE, &lastTimeStep);
    status = H5Aclose(lastAttributeId);

    return lastTimeStep;
}


std::unique_ptr<XFile::TimestepGroup>
XFile::ConcentrationGroup::getTimestepGroup(int timeStep) const {

    std::unique_ptr<XFile::TimestepGroup> tsGroup;
        
    try {
        // Open the sub-group associated with the desired time step.
        tsGroup.reset(new TimestepGroup(*this, timeStep));
    }
    catch(HDF5Exception& e) {
        // We were unable to open the group associated with the given time step.
        assert(not tsGroup);
    }

    return std::move(tsGroup);
}


std::unique_ptr<XFile::TimestepGroup>
XFile::ConcentrationGroup::getLastTimestepGroup(void) const {

    std::unique_ptr<XFile::TimestepGroup> tsGroup;
        
    try {
        // Open the sub-group associated with the last known time step,
        // if any time steps have been written.
        auto lastTimeStep = getLastTimeStep();
        if(lastTimeStep >= 0) {
            tsGroup.reset(new TimestepGroup(*this, lastTimeStep));
        }
    }
    catch(HDF5Exception& e) {
        // We were unable to open the group associated with the given time step.
        assert(not tsGroup);
    }

    return std::move(tsGroup);
}


//----------------------------------------------------------------------------
// TimestepGroup
//
std::string XFile::TimestepGroup::groupNamePrefix = "concentration_";


std::string
XFile::TimestepGroup::makeGroupName(const XFile::ConcentrationGroup& concGroup,
                                        int timeStep) {

	std::ostringstream namestr;
    namestr << concGroup.getName() << '/' << groupNamePrefix << timeStep;
    return namestr.str();
}


XFile::TimestepGroup::TimestepGroup(const XFile::ConcentrationGroup& concGroup,
                                    int timeStep,
                                    double time,
                                    double previousTime,
                                    double deltaTime)
  : HDF5File::Group(concGroup, makeGroupName(concGroup, timeStep), true) {

	// Create, write, and close the absolute time attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(getId(), "absoluteTime",
	H5T_IEEE_F64LE, dataspaceId,
	H5P_DEFAULT, H5P_DEFAULT);
	auto status = H5Awrite(attributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(attributeId);

	// Create, write, and close the previous time attribute
	attributeId = H5Acreate2(getId(), "previousTime", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousTime);
	status = H5Aclose(attributeId);

	// Create, write, and close the timestep time attribute
	attributeId = H5Acreate2(getId(), "deltaTime", H5T_IEEE_F64LE,
			dataspaceId,
			H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(attributeId);
	status = H5Sclose(dataspaceId);
}


XFile::TimestepGroup::TimestepGroup(const XFile::ConcentrationGroup& concGroup,
                                    int timeStep)
  : HDF5File::Group(concGroup, makeGroupName(concGroup, timeStep), false) {

    // Base class opened the group, so nothing else to do.
}


void XFile::TimestepGroup::writeSurface1D(int iSurface,
                                            double nInter,
                                            double previousFlux) {

	// Create, write, and close the surface position attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(getId(), "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	auto status = H5Awrite(attributeId, H5T_STD_I32LE, &iSurface);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of interstitial attribute
	attributeId = H5Acreate2(getId(), "nInterstitial", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nInter);
	status = H5Aclose(attributeId);

	// Create, write, and close the flux of interstitial attribute
	attributeId = H5Acreate2(getId(), "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousFlux);
	status = H5Aclose(attributeId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);

	return;
}


void XFile::TimestepGroup::writeSurface2D(std::vector<int> iSurface,
		std::vector<double> nInter, std::vector<double> previousFlux) {

	// Create the array that will store the indices and fill it
	int size = iSurface.size();
	int indexArray[size];
	for (int i = 0; i < size; i++) {
		indexArray[i] = iSurface[i];
	}

	// Create the dataspace for the dataset with dimension dims
    std::array<hsize_t, 1> dims{ (hsize_t)size };
	hid_t dataspaceId = H5Screate_simple(dims.size(), dims.data(), nullptr);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(getId(), "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write surface array in the dataset
	auto status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &indexArray);

	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the array that will store the quantities and fill it
	double quantityArray[size];
	for (int i = 0; i < size; i++) {
		quantityArray[i] = nInter[i];
	}

	// Create the dataset for the surface indices
	datasetId = H5Dcreate2(getId(), "nInterstitial", H5T_IEEE_F64LE,
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
	datasetId = H5Dcreate2(getId(), "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write quantityArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
}


void XFile::TimestepGroup::writeSurface3D(
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
    std::array<hsize_t, 2> dims{ (hsize_t)xSize, (hsize_t)ySize };
	hid_t dataspaceId = H5Screate_simple(dims.size(), dims.data(), nullptr);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(getId(), "iSurface", H5T_STD_I32LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	auto status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
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
	datasetId = H5Dcreate2(getId(), "nInterstitial", H5T_IEEE_F64LE,
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
	datasetId = H5Dcreate2(getId(), "previousIFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);
	// Close the dataset
	status = H5Dclose(datasetId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);
}


void XFile::TimestepGroup::writeBottom1D(double nHe, 
        double previousHeFlux, double nD,
		double previousDFlux, double nT, double previousTFlux) {

	// Create, write, and close the quantity of helium attribute
	hid_t dataspaceId = H5Screate(H5S_SCALAR);
	hid_t attributeId = H5Acreate2(getId(), "nHelium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	auto status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nHe);
	status = H5Aclose(attributeId);
	// Create, write, and close the flux of helium attribute
	attributeId = H5Acreate2(getId(), "previousHeFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousHeFlux);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of deuterium attribute
	attributeId = H5Acreate2(getId(), "nDeuterium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nD);
	status = H5Aclose(attributeId);
	// Create, write, and close the flux of deuterium attribute
	attributeId = H5Acreate2(getId(), "previousDFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousDFlux);
	status = H5Aclose(attributeId);

	// Create, write, and close the quantity of tritium attribute
	attributeId = H5Acreate2(getId(), "nTritium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &nT);
	status = H5Aclose(attributeId);
	// Create, write, and close the flux of tritium attribute
	attributeId = H5Acreate2(getId(), "previousTFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attributeId, H5T_IEEE_F64LE, &previousTFlux);
	status = H5Aclose(attributeId);

	// Close the dataspace
	status = H5Sclose(dataspaceId);
}


void XFile::TimestepGroup::writeBottom2D(std::vector<double> nHe,
		std::vector<double> previousHeFlux, std::vector<double> nD,
		std::vector<double> previousDFlux, std::vector<double> nT,
		std::vector<double> previousTFlux) {

	// Find out the size of the arrays
	const int size = nHe.size();

	// Create the dataspace for the dataset with dimension dims
    std::array<hsize_t, 1> dims{ (hsize_t)size };
	hid_t dataspaceId = H5Screate_simple(dims.size(), dims.data(), nullptr);

	// Create the dataset for helium
	hid_t datasetId = H5Dcreate2(getId(), "nHelium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	// Create the array that will store the quantities and fill it
//	double quantityArray[size];
//	for (int i = 0; i < size; i++) {
//		quantityArray[i] = nInter[i];
//	}
	// Write the dataset
	auto status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nHe.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the helium flux
	datasetId = H5Dcreate2(getId(), "previousHeFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousHeFlux.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the deuterium
	datasetId = H5Dcreate2(getId(), "nDeuterium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nD.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the deuterium flux
	datasetId = H5Dcreate2(getId(), "previousDFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousDFlux.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the tritium
	datasetId = H5Dcreate2(getId(), "nTritium", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nT.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the tritium flux
	datasetId = H5Dcreate2(getId(), "previousTFlux", H5T_IEEE_F64LE,
			dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousTFlux.data());

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);
}


void XFile::TimestepGroup::writeConcentrationDataset(int size,
        double concArray[][2],
        int i, int j, int k) {

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Create the dataspace for the dataset with dimension dims
    std::array<hsize_t, 2> dims{ (hsize_t)size, (hsize_t)2 };
	auto concDataspaceId = H5Screate_simple(dims.size(), dims.data(), nullptr);

	// Create the dataset of concentrations for this position
	hid_t datasetId = H5Dcreate2(getId(), datasetName.str().c_str(),
	H5T_IEEE_F64LE, concDataspaceId,
	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Close dataspace
	auto status = H5Sclose(concDataspaceId);

	// Write concArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, concArray);

	// Close dataset
	status = H5Dclose(datasetId);

	return;
}


std::pair<double, double> XFile::TimestepGroup::readTimes(void) const {

    double time;
    double deltaTime;

	// Open and read the absoluteTime attribute
	hid_t attributeId = H5Aopen(getId(), "absoluteTime", H5P_DEFAULT);
	auto status = H5Aread(attributeId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(attributeId);

	// Open and read the deltaTime attribute
	attributeId = H5Aopen(getId(), "deltaTime", H5P_DEFAULT);
	status = H5Aread(attributeId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(attributeId);

    return std::make_pair(time, deltaTime);
}


double XFile::TimestepGroup::readPreviousTime(void) const {

	// Open and read the previousTime attribute
	double previousTime = 0.0;
	hid_t attributeId = H5Aopen(getId(), "previousTime", H5P_DEFAULT);
	auto status = H5Aread(attributeId, H5T_IEEE_F64LE, &previousTime);

	// Close everything
	status = H5Aclose(attributeId);

	return previousTime;
}


int XFile::TimestepGroup::readSurface1D(void) const {

	// Initialize the surface position
	int iSurface = 0;

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(getId(), "iSurface", H5P_DEFAULT);
	auto status = H5Aread(attributeId, H5T_STD_I32LE, &iSurface);
	status = H5Aclose(attributeId);

	return iSurface;
}

std::vector<int> XFile::TimestepGroup::readSurface2D(void) const {

	// Create the vector to return
	std::vector<int> toReturn;

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), "iSurface", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
    std::array<hsize_t, 1> dims;
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

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

	return toReturn;
}

std::vector<std::vector<int> > XFile::TimestepGroup::readSurface3D(void) const {

	// Create the vector to return
	std::vector<std::vector<int> > toReturn;

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), "iSurface", H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
    std::array<hsize_t, 2> dims;
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

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

	return toReturn;
}


double XFile::TimestepGroup::readData1D(const std::string& dataName) const {

	// Initialize the surface position
	double data = 0.0;

	// Open and read the iSurface attribute
	hid_t attributeId = H5Aopen(getId(), dataName.c_str(), H5P_DEFAULT);
	auto status = H5Aread(attributeId, H5T_IEEE_F64LE, &data);
	status = H5Aclose(attributeId);

	return data;
}

std::vector<double> XFile::TimestepGroup::readData2D(
                                const std::string& dataName) const {

	// Create the vector to return
	std::vector<double> toReturn;

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), dataName.c_str(), H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
    std::array<hsize_t, 1> dims;
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

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

	return toReturn;
}

std::vector<std::vector<double> > XFile::TimestepGroup::readData3D(
		const std::string& dataName) const {

	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), dataName.c_str(), H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
    std::array<hsize_t, 2> dims;
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

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

	return toReturn;
}


std::vector<std::vector<double> > XFile::TimestepGroup::readGridPoint(
        int i, int j, int k) const {

	// Create the vector to return
	std::vector<std::vector<double> > toReturn;

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Check the dataset
	bool datasetExist = H5Lexists(getId(), datasetName.str().c_str(),
	H5P_DEFAULT);

	// If the dataset exists
	if (datasetExist) {
		// Open the dataset
		hid_t datasetId = H5Dopen(getId(), datasetName.str().c_str(),
		H5P_DEFAULT);

		// Get the dataspace object
		hid_t dataspaceId = H5Dget_space(datasetId);

		// Get the dimensions of the dataset
        std::array<hsize_t, 2> dims;
		auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

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

	return toReturn;
}


} // namespace xolotlCore

