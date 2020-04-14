#include <iostream>
#include <sstream>
#include <iterator>
#include <array>
#include "hdf5.h"
#include "mpi.h"
#include "xolotlCore/io/XFile.h"
#include <PSISuperCluster.h>
#include <FeSuperCluster.h>
#include <NESuperCluster.h>
#include <AlloySuperCluster.h>

namespace xolotlCore {

HDF5File::AccessMode XFile::EnsureCreateAccessMode(HDF5File::AccessMode mode) {

	bool createMode = ((mode == HDF5File::AccessMode::CreateOrTruncateIfExists)
			or (mode == HDF5File::AccessMode::CreateOrFailIfExists));

	if (not createMode) {
		throw HDF5Exception(
				"Attempt to create existing file with a non-create access mode");
	}
	return mode;
}

HDF5File::AccessMode XFile::EnsureOpenAccessMode(HDF5File::AccessMode mode) {

	bool openMode = ((mode == HDF5File::AccessMode::OpenReadOnly)
			or (mode == HDF5File::AccessMode::OpenReadWrite));

	if (not openMode) {
		throw HDF5Exception(
				"Attempt to open existing file with a non-open access mode");
	}
	return mode;
}

XFile::XFile(fs::path _path, const std::vector<double> &grid, MPI_Comm _comm,
		int ny, double hy, int nz, double hz, AccessMode _mode) :
		HDF5File(_path, EnsureCreateAccessMode(_mode), _comm, true) {
	// Create and initialize the header group.
	HeaderGroup headerGroup(*this, grid, ny, hy, nz, hz);

	// Create and initialize the group where the concentrations will be stored
	ConcentrationGroup concGroup(*this, true);
}

XFile::XFile(fs::path _path, MPI_Comm _comm, AccessMode _mode) :
		HDF5File(_path, EnsureOpenAccessMode(_mode), _comm, true) {

	// Nothing else to do.
}

//----------------------------------------------------------------------------
// HeaderGroup
//
const fs::path XFile::HeaderGroup::path = "/headerGroup";
const std::string XFile::HeaderGroup::nxAttrName = "nx";
const std::string XFile::HeaderGroup::hxAttrName = "hx";
const std::string XFile::HeaderGroup::nyAttrName = "ny";
const std::string XFile::HeaderGroup::hyAttrName = "hy";
const std::string XFile::HeaderGroup::nzAttrName = "nz";
const std::string XFile::HeaderGroup::hzAttrName = "hz";

XFile::HeaderGroup::HeaderGroup(const XFile &file,
		const std::vector<double> &grid, int ny, double hy, int nz, double hz) :
		HDF5File::Group(file, HeaderGroup::path, true) {

	// Base class created the group.

	// Build a dataspace for our scalar attributes.
	XFile::ScalarDataSpace scalarDSpace;

	// Add an nx attribute.
	int nx = grid.size() - 2;
	Attribute<decltype(nx)> nxAttr(*this, nxAttrName, scalarDSpace);
	nxAttr.setTo(nx);

	// Add an hy attribute.
	double hx = 0.0;
	if (grid.size() > 0)
		hx = grid[1] - grid[0];
	Attribute<decltype(hx)> hxAttr(*this, hxAttrName, scalarDSpace);
	hxAttr.setTo(hx);

	// Add an ny attribute.
	Attribute<decltype(ny)> nyAttr(*this, nyAttrName, scalarDSpace);
	nyAttr.setTo(ny);
	// Add an hy attribute.
	Attribute<decltype(hy)> hyAttr(*this, hyAttrName, scalarDSpace);
	hyAttr.setTo(hy);

	// Add an nz attribute.
	Attribute<decltype(nz)> nzAttr(*this, nzAttrName, scalarDSpace);
	nzAttr.setTo(nz);
	// Add an hz attribute.
	Attribute<decltype(hz)> hzAttr(*this, hzAttrName, scalarDSpace);
	hzAttr.setTo(hz);

	if (nx > 0) {
		// Create, write, and close the grid dataset
		double gridArray[nx];
		for (int i = 0; i < nx; i++) {
			gridArray[i] = grid[i + 1] - grid[1];
		}
		std::array<hsize_t, 1> dims { (hsize_t) nx };
		XFile::SimpleDataSpace<1> gridDSpace(dims);
		hid_t datasetId = H5Dcreate2(getId(), "grid", H5T_IEEE_F64LE,
				gridDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		auto status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, &gridArray);
		status = H5Dclose(datasetId);
	}
}

XFile::HeaderGroup::HeaderGroup(const XFile &file) :
		HDF5File::Group(file, HeaderGroup::path, false) {

	// Base class opened the group, so nothing else to do.
}

void XFile::HeaderGroup::read(int &nx, double &hx, int &ny, double &hy, int &nz,
		double &hz) const {

	// Open and read the nx attribute
	Attribute<int> nxAttr(*this, nxAttrName);
	nx = nxAttr.get();

	// Open and read the hx attribute
	Attribute<double> hxAttr(*this, hxAttrName);
	hx = hxAttr.get();

	// Open and read the ny attribute
	Attribute<int> nyAttr(*this, nyAttrName);
	ny = nyAttr.get();

	// Open and read the hy attribute
	Attribute<double> hyAttr(*this, hyAttrName);
	hy = hyAttr.get();

	// Open and read the nz attribute
	Attribute<int> nzAttr(*this, nzAttrName);
	nz = nzAttr.get();

	// Open and read the hz attribute
	Attribute<double> hzAttr(*this, hzAttrName);
	hz = hzAttr.get();
}

//----------------------------------------------------------------------------
// NetworkGroup
//
const fs::path XFile::NetworkGroup::path = "/networkGroup";
const std::string XFile::NetworkGroup::sizeAttrName = "totalSize";
const std::string XFile::NetworkGroup::phaseSpaceAttrName = "phaseSpace";

XFile::NetworkGroup::NetworkGroup(const XFile &file) :
		HDF5File::Group(file, NetworkGroup::path, false) {

	// Base class opened the group, so nothing else to do.
}

XFile::NetworkGroup::NetworkGroup(const XFile &file,
		experimental::IReactionNetwork &network) :
		HDF5File::Group(file, NetworkGroup::path, true) {
	// Base class created the group.

	// Get the phase space information
	auto phaseSpace = network.getPhaseSpace();
	// Write it as an attribute
//	H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
//	std::array<hsize_t, 1> dim { phaseSpace.size() };
//	XFile::SimpleDataSpace<1> phaseDSpace(dim);
//	hid_t attrId = H5Acreate2(getId(), phaseSpaceAttrName.c_str(), datatype,
//			phaseDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT);
//	auto status = H5Awrite(attrId, datatype, phaseSpace.data());
//	status = H5Aclose(attrId);

	// Get the bounds for each cluster
	auto bounds = network.getAllClusterBounds();
	// Get the size information
	int totalSize = network.getNumClusters();

	// Build a dataspace for our scalar attributes.
	XFile::ScalarDataSpace scalarDSpace;

	// Add a total size attribute.
	Attribute<decltype(totalSize)> normalSizeAttr(*this, sizeAttrName,
			scalarDSpace);
	normalSizeAttr.setTo(totalSize);

	for (experimental::IReactionNetwork::IndexType i = 0; i < totalSize; i++) {
		auto cluster = network.getClusterCommon(i);
		// Create and initialize the cluster group
		ClusterGroup clusterGroup(*this, i, bounds[i],
				cluster.getFormationEnergy(), cluster.getMigrationEnergy(),
				cluster.getDiffusionFactor());
	}
}

int XFile::NetworkGroup::readNetworkSize() const {
	// Open and read the total size attribute
	Attribute<int> totalSizeAttr(*this, sizeAttrName);
	return totalSizeAttr.get();
}

void XFile::NetworkGroup::readReactions(
		experimental::IReactionNetwork &network) const {
	// Doesn't do anything for now

	return;
}

void XFile::NetworkGroup::copyTo(const XFile &target) const {

	H5Ocopy(getLocation().getId(), NetworkGroup::path.string().c_str(),
			target.getId(), NetworkGroup::path.string().c_str(),
			H5P_DEFAULT,
			H5P_DEFAULT);
}

//----------------------------------------------------------------------------
// ClusterGroup
//
const std::string XFile::ClusterGroup::formationEnergyAttrName =
		"formationEnergy";
const std::string XFile::ClusterGroup::migrationEnergyAttrName =
		"migrationEnergy";
const std::string XFile::ClusterGroup::diffusionFactorAttrName =
		"diffusionFactor";
const std::string XFile::ClusterGroup::boundsAttrName = "bounds";

XFile::ClusterGroup::ClusterGroup(const NetworkGroup &networkGroup, int id) :
		HDF5File::Group(networkGroup, makeGroupName(id), false) {
}

XFile::ClusterGroup::ClusterGroup(const NetworkGroup &networkGroup, int id,
		ClusterBoundsType bounds, double formationEnergy,
		double migrationEnergy, double diffusionFactor) :
		HDF5File::Group(networkGroup, makeGroupName(id), true) {
	// Write the region
	std::array<hsize_t, 1> dim { bounds.size() };
	XFile::SimpleDataSpace<1> boundDSpace(dim);
	hid_t attrId = H5Acreate2(getId(), boundsAttrName.c_str(),
	H5T_STD_I32LE, boundDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT);
	auto status = H5Awrite(attrId, H5T_STD_I32LE, bounds.data());
	status = H5Aclose(attrId);

	// Build a dataspace for our scalar attributes.
	XFile::ScalarDataSpace scalarDSpace;

	// Add a formationEnergy attribute.
	Attribute<decltype(formationEnergy)> formationEnergyAttr(*this,
			formationEnergyAttrName, scalarDSpace);
	formationEnergyAttr.setTo(formationEnergy);
	// Add a migrationEnergy attribute.
	Attribute<decltype(migrationEnergy)> migrationEnergyAttr(*this,
			migrationEnergyAttrName, scalarDSpace);
	migrationEnergyAttr.setTo(migrationEnergy);
	// Add a diffusionFactor attribute.
	Attribute<decltype(diffusionFactor)> diffusionFactorAttr(*this,
			diffusionFactorAttrName, scalarDSpace);
	diffusionFactorAttr.setTo(diffusionFactor);

	return;
}

std::string XFile::ClusterGroup::makeGroupName(int id) {
	std::ostringstream namestr;
	namestr << id;
	return namestr.str();
}

XFile::ClusterGroup::ClusterBoundsType XFile::ClusterGroup::readCluster(
		double &formationEnergy, double &migrationEnergy,
		double &diffusionFactor) const {
	// Open and read the formation energy attribute
	Attribute<double> formationEnergyAttr(*this, formationEnergyAttrName);
	formationEnergy = formationEnergyAttr.get();
	// Open and read the migration energy attribute
	Attribute<double> migrationEnergyAttr(*this, migrationEnergyAttrName);
	migrationEnergy = migrationEnergyAttr.get();
	// Open and read the diffusion factor attribute
	Attribute<double> diffusionFactorAttr(*this, diffusionFactorAttrName);
	diffusionFactor = diffusionFactorAttr.get();

	// Read the composition attribute
	hid_t attributeId = H5Aopen_name(getId(), boundsAttrName.c_str());
	hid_t dataspaceId = H5Aget_space(attributeId);
	hsize_t dims[1];
	herr_t status = H5Sget_simple_extent_dims(dataspaceId, dims, NULL);
	ClusterBoundsType bounds;
	bounds.resize(dims[0]);
	status = H5Aread(attributeId, H5T_STD_I32LE, bounds.data());
	status = H5Aclose(attributeId);

	return bounds;
}

//----------------------------------------------------------------------------
// ConcentrationGroup
//
const fs::path XFile::ConcentrationGroup::path = "/concentrationsGroup";
const std::string XFile::ConcentrationGroup::lastTimestepAttrName =
		"lastTimeStep";

XFile::ConcentrationGroup::ConcentrationGroup(const XFile &file, bool create) :
		HDF5File::Group(file, ConcentrationGroup::path, create) {

	if (create) {
		// We just created the group.
		// Populate it.
		//
		// Create, write, and close the last written time step attribute
		int lastTimeStep = -1;
		XFile::ScalarDataSpace lastDSpace;
		Attribute<decltype(lastTimeStep)> lastTimestepAttr(*this,
				lastTimestepAttrName, lastDSpace);
		lastTimestepAttr.setTo(lastTimeStep);
	}
}

std::unique_ptr<XFile::TimestepGroup> XFile::ConcentrationGroup::addTimestepGroup(
		int timeStep, double time, double previousTime,
		double deltaTime) const {

	// Create a group for the new timestep.
	std::unique_ptr<XFile::TimestepGroup> tsGroup(
			new TimestepGroup(*this, timeStep, time, previousTime, deltaTime));

	// Update our last known timestep.
	Attribute<decltype(timeStep)> lastTimestepAttr(*this, lastTimestepAttrName);
	lastTimestepAttr.setTo(timeStep);

	return std::move(tsGroup);
}

int XFile::ConcentrationGroup::getLastTimeStep(void) const {

	Attribute<int> lastTimestepAttr(*this, lastTimestepAttrName);
	return lastTimestepAttr.get();
}

std::unique_ptr<XFile::TimestepGroup> XFile::ConcentrationGroup::getTimestepGroup(
		int timeStep) const {

	std::unique_ptr<XFile::TimestepGroup> tsGroup;

	try {
		// Open the sub-group associated with the desired time step.
		tsGroup.reset(new TimestepGroup(*this, timeStep));
	} catch (HDF5Exception &e) {
		// We were unable to open the group associated with the given time step.
		assert(not tsGroup);
	}

	return std::move(tsGroup);
}

std::unique_ptr<XFile::TimestepGroup> XFile::ConcentrationGroup::getLastTimestepGroup(
		void) const {

	std::unique_ptr<XFile::TimestepGroup> tsGroup;

	try {
		// Open the sub-group associated with the last known time step,
		// if any time steps have been written.
		auto lastTimeStep = getLastTimeStep();
		if (lastTimeStep >= 0) {
			tsGroup.reset(new TimestepGroup(*this, lastTimeStep));
		}
	} catch (HDF5Exception &e) {
		// We were unable to open the group associated with the given time step.
		assert(not tsGroup);
	}

	return std::move(tsGroup);
}

//----------------------------------------------------------------------------
// TimestepGroup
//
const std::string XFile::TimestepGroup::groupNamePrefix = "concentration_";
const std::string XFile::TimestepGroup::absTimeAttrName = "absoluteTime";
const std::string XFile::TimestepGroup::prevTimeAttrName = "previousTime";
const std::string XFile::TimestepGroup::deltaTimeAttrName = "deltaTime";
const std::string XFile::TimestepGroup::surfacePosDataName = "iSurface";
const std::string XFile::TimestepGroup::nIntersAttrName = "nInterstitial";
const std::string XFile::TimestepGroup::prevIFluxAttrName = "previousIFlux";
const std::string XFile::TimestepGroup::nHeAttrName = "nHelium";
const std::string XFile::TimestepGroup::prevHeFluxAttrName = "previousHeFlux";
const std::string XFile::TimestepGroup::nDAttrName = "nDeuterium";
const std::string XFile::TimestepGroup::prevDFluxAttrName = "previousDFlux";
const std::string XFile::TimestepGroup::nTAttrName = "nTritium";
const std::string XFile::TimestepGroup::prevTFluxAttrName = "previousTFlux";
const std::string XFile::TimestepGroup::nVAttrName = "nVacancy";
const std::string XFile::TimestepGroup::prevVFluxAttrName = "previousVFlux";
const std::string XFile::TimestepGroup::nIntersBulkAttrName = "nIBulk";
const std::string XFile::TimestepGroup::prevIBulkFluxAttrName =
		"previousIBulkFlux";

const std::string XFile::TimestepGroup::concDatasetName = "concs";

std::string XFile::TimestepGroup::makeGroupName(
		const XFile::ConcentrationGroup &concGroup, int timeStep) {

	std::ostringstream namestr;
	namestr << concGroup.getName() << '/' << groupNamePrefix << timeStep;
	return namestr.str();
}

XFile::TimestepGroup::TimestepGroup(const XFile::ConcentrationGroup &concGroup,
		int timeStep, double time, double previousTime, double deltaTime) :
		HDF5File::Group(concGroup, makeGroupName(concGroup, timeStep), true) {

	// Get a dataspace for our scalar attributes.
	XFile::ScalarDataSpace scalarDSpace;

	// Add absolute time attribute.
	Attribute<decltype(time)> absTimeAttr(*this, absTimeAttrName, scalarDSpace);
	absTimeAttr.setTo(time);

	// Add previous time attribute.
	Attribute<decltype(previousTime)> prevTimeAttr(*this, prevTimeAttrName,
			scalarDSpace);
	prevTimeAttr.setTo(previousTime);

	// Add delta time attribute.
	Attribute<decltype(deltaTime)> deltaTimeAttr(*this, deltaTimeAttrName,
			scalarDSpace);
	deltaTimeAttr.setTo(deltaTime);
}

XFile::TimestepGroup::TimestepGroup(const XFile::ConcentrationGroup &concGroup,
		int timeStep) :
		HDF5File::Group(concGroup, makeGroupName(concGroup, timeStep), false) {

	// Base class opened the group, so nothing else to do.
}

void XFile::TimestepGroup::writeSurface1D(Surface1DType iSurface,
		Data1DType nInter, Data1DType previousFlux) const {

	// Make a scalar dataspace for 1D attributes.
	XFile::ScalarDataSpace scalarDSpace;

	// Create, write, and close the surface position attribute
	Attribute<int> surfacePosAttr(*this, surfacePosDataName, scalarDSpace);
	surfacePosAttr.setTo(iSurface);

	// Create, write, and close the quantity of interstitial attribute
	Attribute<double> nIntersAttr(*this, nIntersAttrName, scalarDSpace);
	nIntersAttr.setTo(nInter);

	// Create, write, and close the flux of interstitial attribute
	Attribute<double> prevIFluxAttr(*this, prevIFluxAttrName, scalarDSpace);
	prevIFluxAttr.setTo(previousFlux);

	return;
}

void XFile::TimestepGroup::writeSurface2D(const Surface2DType &iSurface,
		const Data2DType &nInter, const Data2DType &previousFlux) const {

	// Create the array that will store the indices and fill it
	int size = iSurface.size();

	// Create the dataspace for the dataset with dimension dims
	std::array<hsize_t, 1> dims { (hsize_t) size };
	XFile::SimpleDataSpace<1> indexDSpace(dims);

	// Create and write the dataset for the surface indices
	auto datasetId = H5Dcreate2(getId(), surfacePosDataName.c_str(),
	H5T_STD_I32LE, indexDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);

	// Write surface index data in the dataset
	auto status = H5Dwrite(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, iSurface.data());

	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the array that will store the quantities and fill it
	double quantityArray[size];
	for (int i = 0; i < size; i++) {
		quantityArray[i] = nInter[i];
	}

	// Create the dataset for the surface indices
	datasetId = H5Dcreate2(getId(), nIntersAttrName.c_str(), H5T_IEEE_F64LE,
			indexDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
	datasetId = H5Dcreate2(getId(), prevIFluxAttrName.c_str(), H5T_IEEE_F64LE,
			indexDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write quantityArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);

	// Close everything
	status = H5Dclose(datasetId);
}

void XFile::TimestepGroup::writeSurface3D(const Surface3DType &iSurface,
		const Data3DType &nInter, const Data3DType &previousFlux) const {

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
	std::array<hsize_t, 2> dims { (hsize_t) xSize, (hsize_t) ySize };
	XFile::SimpleDataSpace<2> indexDSpace(dims);

	// Create the dataset for the surface indices
	hid_t datasetId = H5Dcreate2(getId(), surfacePosDataName.c_str(),
	H5T_STD_I32LE, indexDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT,
	H5P_DEFAULT);
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
	datasetId = H5Dcreate2(getId(), nIntersAttrName.c_str(), H5T_IEEE_F64LE,
			indexDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
	datasetId = H5Dcreate2(getId(), prevIFluxAttrName.c_str(), H5T_IEEE_F64LE,
			indexDSpace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &quantityArray);
	// Close the dataset
	status = H5Dclose(datasetId);
}

void XFile::TimestepGroup::writeBottom1D(Data1DType nHe,
		Data1DType previousHeFlux, Data1DType nD, Data1DType previousDFlux,
		Data1DType nT, Data1DType previousTFlux, Data1DType nV,
		Data1DType previousVFlux, Data1DType nI, Data1DType previousIFlux) {

	// Build a data space for scalar attributes.
	XFile::ScalarDataSpace scalarDSpace;

	// Add quantity of helium attribute
	Attribute<Data1DType> nHeAttr(*this, nHeAttrName, scalarDSpace);
	nHeAttr.setTo(nHe);

	// Add flux of helium attribute
	Attribute<Data1DType> prevHeFluxAttr(*this, prevHeFluxAttrName,
			scalarDSpace);
	prevHeFluxAttr.setTo(previousHeFlux);

	// Add quantity of deuterium attribute
	Attribute<Data1DType> nDAttr(*this, nDAttrName, scalarDSpace);
	nDAttr.setTo(nD);

	// Add flux of deuterium attribute
	Attribute<Data1DType> prevDFluxAttr(*this, prevDFluxAttrName, scalarDSpace);
	prevDFluxAttr.setTo(previousDFlux);

	// Add quantity of tritium attribute
	Attribute<Data1DType> nTAttr(*this, nTAttrName, scalarDSpace);
	nTAttr.setTo(nT);

	// Add flux of tritium attribute
	Attribute<Data1DType> prevTFluxAttr(*this, prevTFluxAttrName, scalarDSpace);
	prevTFluxAttr.setTo(previousTFlux);

	// Add quantity of vacancy attribute
	Attribute<Data1DType> nVAttr(*this, nVAttrName, scalarDSpace);
	nVAttr.setTo(nV);

	// Add flux of vacancy attribute
	Attribute<Data1DType> prevVFluxAttr(*this, prevVFluxAttrName, scalarDSpace);
	prevVFluxAttr.setTo(previousVFlux);

	// Add quantity of int attribute
	Attribute<Data1DType> nIAttr(*this, nIntersBulkAttrName, scalarDSpace);
	nIAttr.setTo(nI);

	// Add flux of int attribute
	Attribute<Data1DType> prevIFluxAttr(*this, prevIBulkFluxAttrName,
			scalarDSpace);
	prevIFluxAttr.setTo(previousIFlux);
}

void XFile::TimestepGroup::writeBottom2D(const Data2DType &nHe,
		const Data2DType &previousHeFlux, const Data2DType &nD,
		const Data2DType &previousDFlux, const Data2DType &nT,
		const Data2DType &previousTFlux) {

	// Find out the size of the arrays
	const int size = nHe.size();

	// Create the dataspace for the dataset with dimension dims
	std::array<hsize_t, 1> dims { (hsize_t) size };
	XFile::SimpleDataSpace<1> dspace(dims);

	// Create the dataset for helium
	hid_t datasetId = H5Dcreate2(getId(), "nHelium", H5T_IEEE_F64LE,
			dspace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
			dspace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousHeFlux.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the deuterium
	datasetId = H5Dcreate2(getId(), "nDeuterium", H5T_IEEE_F64LE,
			dspace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nD.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the deuterium flux
	datasetId = H5Dcreate2(getId(), "previousDFlux", H5T_IEEE_F64LE,
			dspace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousDFlux.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the tritium
	datasetId = H5Dcreate2(getId(), "nTritium", H5T_IEEE_F64LE, dspace.getId(),
	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, nT.data());
	// Close the dataset
	status = H5Dclose(datasetId);

	// Create the dataset for the tritium flux
	datasetId = H5Dcreate2(getId(), "previousTFlux", H5T_IEEE_F64LE,
			dspace.getId(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write  the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, previousTFlux.data());

	// Close everything
	status = H5Dclose(datasetId);
}

void XFile::TimestepGroup::writeConcentrationDataset(int size,
		double concArray[][2], bool write, int i, int j, int k) {

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Create the dataspace for the dataset with dimension dims
	std::array<hsize_t, 2> dims { (hsize_t) size, (hsize_t) 2 };
	XFile::SimpleDataSpace<2> concDSpace(dims);

	// Create the dataset of concentrations for this position
	hid_t datasetId = H5Dcreate2(getId(), datasetName.str().c_str(),
	H5T_IEEE_F64LE, concDSpace.getId(),
	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create property list for independent dataset write.
	hid_t propertyListId = H5Pcreate(H5P_DATASET_XFER);
	auto status = H5Pset_dxpl_mpio(propertyListId, H5FD_MPIO_INDEPENDENT);

	if (write) {
		// Write concArray in the dataset
		status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
				propertyListId, concArray);
	}

	// Close dataset
	status = H5Dclose(datasetId);

	return;
}

// Caller gives us 2D ragged representation, and we flatten it into
// a 1D dataset and add a 1D "starting index" array.
// Assumes that grid point slabs are assigned to processes in 
// MPI rank order.
void XFile::TimestepGroup::writeConcentrations(const XFile &file, int baseX,
		const Concs1DType &raggedConcs) const {

	// Create and write the ragged dataset.
	RaggedDataSet2D<ConcType> dataset(file.getComm(), *this, concDatasetName,
			baseX, raggedConcs);

	// Unlike our other DataSet types, there is no need to call a
	// 'write' on the dataset.  The constructor above
	// defines the dataset *and* writes the given data.
}

XFile::TimestepGroup::Concs1DType XFile::TimestepGroup::readConcentrations(
		const XFile &file, int baseX, int numX) const {

	// Open and read the ragged dataset.
	RaggedDataSet2D<ConcType> dataset(file.getComm(), *this, concDatasetName);
	return dataset.read(baseX, numX);
}

std::pair<double, double> XFile::TimestepGroup::readTimes(void) const {

	// Open the desired attributes.
	Attribute<double> absTimeAttr(*this, absTimeAttrName);
	Attribute<double> deltaTimeAttr(*this, deltaTimeAttrName);

	// Read and return the attributes.
	return std::make_pair(absTimeAttr.get(), deltaTimeAttr.get());
}

double XFile::TimestepGroup::readPreviousTime(void) const {

	// Open and read the previousTime attribute
	Attribute<double> prevTimeAttr(*this, prevTimeAttrName);
	return prevTimeAttr.get();
}

auto XFile::TimestepGroup::readSurface1D(void) const -> Surface1DType {

	// Open and read the surface position attribute.
	Attribute<Surface1DType> surfacePosAttr(*this, surfacePosDataName);
	return surfacePosAttr.get();
}

auto XFile::TimestepGroup::readSurface2D(void) const -> Surface2DType {

	DataSet<Surface2DType> dataset(*this, surfacePosDataName);
	return dataset.read();
}

auto XFile::TimestepGroup::readSurface3D(void) const -> Surface3DType {

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), surfacePosDataName.c_str(), H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	std::array<hsize_t, 2> dims;
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

	// Create the array that will receive the indices
	Surface3DType::value_type index(dims[0] * dims[1]);

	// Read the data set
	status = H5Dread(datasetId, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, index.data());

	// Loop on the length and fill the vector to return
	Surface3DType toReturn(dims[0]);
	for (int i = 0; i < dims[0]; ++i) {

		toReturn[i].reserve(dims[1]);
		std::copy(index.begin() + i * dims[1],
				index.begin() + (i + 1) * dims[1],
				std::back_inserter(toReturn[i]));
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);

	return toReturn;
}

auto XFile::TimestepGroup::readData1D(const std::string &dataName) const ->
		Data1DType {

	Attribute<Data1DType> attr(*this, dataName);
	return attr.get();

}

auto XFile::TimestepGroup::readData2D(const std::string &dataName) const ->
		Data2DType {

	DataSet<Data2DType> dataset(*this, dataName);
	return dataset.read();
}

auto XFile::TimestepGroup::readData3D(const std::string &dataName) const ->
		Data3DType {

	// Open the dataset
	hid_t datasetId = H5Dopen(getId(), dataName.c_str(), H5P_DEFAULT);

	// Get the dataspace object
	hid_t dataspaceId = H5Dget_space(datasetId);

	// Get the dimensions of the dataset
	std::array<hsize_t, 2> dims;
	auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(), nullptr);

	// Create the array that will receive the indices
	Data3DType::value_type quantity(dims[0] * dims[1]);

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, quantity.data());

	// Loop on the length and fill the vector to return
	Data3DType toReturn(dims[0]);
	for (int i = 0; i < dims[0]; i++) {

		toReturn[i].reserve(dims[1]);
		std::copy(quantity.begin() + i * dims[1],
				quantity.begin() + (i + 1) * dims[1],
				std::back_inserter(toReturn[i]));
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Sclose(dataspaceId);

	return toReturn;
}

auto XFile::TimestepGroup::readGridPoint(int i, int j, int k) const ->
		Data3DType {

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << i << "_" << j << "_" << k;

	// Check the dataset
	bool datasetExist = H5Lexists(getId(), datasetName.str().c_str(),
	H5P_DEFAULT);

	// If the dataset exists
	Data3DType toReturn;
	if (datasetExist) {
		// Open the dataset
		hid_t datasetId = H5Dopen(getId(), datasetName.str().c_str(),
		H5P_DEFAULT);

		// Get the dataspace object
		hid_t dataspaceId = H5Dget_space(datasetId);

		// Get the dimensions of the dataset
		std::array<hsize_t, 2> dims;
		auto status = H5Sget_simple_extent_dims(dataspaceId, dims.data(),
				nullptr);

		// Create the array that will receive the concentrations
		Data3DType::value_type conc(dims[0] * dims[1]);

		// Read the data set
		status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, conc.data());

		// Loop on the length
		toReturn.resize(dims[0]);
		for (unsigned int n = 0; n < dims[0]; n++) {

			toReturn[n].reserve(dims[1]);
			std::copy(conc.begin() + n * dims[1],
					conc.begin() + (n + 1) * dims[1],
					std::back_inserter(toReturn[n]));
		}

		// Close everything
		status = H5Dclose(datasetId);
		status = H5Sclose(dataspaceId);
	}

	return toReturn;
}

} // namespace xolotlCore

