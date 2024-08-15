#ifndef XCORE_XFILE_H
#define XCORE_XFILE_H

#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <xolotl/core/network/IReactionNetwork.h>
#include <xolotl/io/HDF5Exception.h>
#include <xolotl/io/HDF5File.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace io
{
// Class for reading and writing an HDF5 file with Xolotl data.
// Note: the class stores 1D data as an attribute on a group instead
// of as a dataset.
// TODO Why?  because it is an attribute, every process must have same
// data to write.  As opposed to dataset, where can use independent file
// access to let one process write. (?)
class XFile : public HDF5File
{
public:
	using NetworkType = core::network::IReactionNetwork;

	// A group with info about a specific time step.
	class ConcentrationGroup;
	class TimestepGroup : public HDF5File::Group
	{
	private:
		// Prefix to use when constructing group names.
		static const std::string groupNamePrefix;

		// Names of time-related attributes.
		static const std::string absTimeAttrName;
		static const std::string prevTimeAttrName;
		static const std::string deltaTimeAttrName;

		// Names of surface position attributes.
		static const std::string surfacePosDataName;

		// Names of Helium attributes.
		static const std::string nHeBurstAttrName;

		// Names of Deuterium attributes.
		static const std::string nDBurstAttrName;

		// Names of Tritium attributes.
		static const std::string nTBurstAttrName;

		// Names of surface and bulk attributes.
		static const std::string nAttrName;
		static const std::string previousFluxAttrName;
		static const std::string surfAttrName;
		static const std::string bulkAttrName;

		// Name of the concentrations data set.
		static const std::string concDatasetName;

		// Names of grid-specification attributes.
		static const std::string nxAttrName;
		static const std::string hxAttrName;
		static const std::string nyAttrName;
		static const std::string hyAttrName;
		static const std::string nzAttrName;
		static const std::string hzAttrName;

	public:
		// Concise name for surface representations.
		using Surface1DType = int;
		using Surface2DType = std::vector<Surface1DType>;
		using Surface3DType = std::vector<Surface2DType>;

		// Concise name for data representations.
		using Data1DType = double;
		using Data2DType = std::vector<Data1DType>;
		using Data3DType = std::vector<Data2DType>;

		// Concise name for concentrations data type.
		// Because the number of concentrations we write for each
		// grid point can vary, these multidimensional data types must
		// support ragged edges in the last dimension.
		using ConcType = std::pair<int, double>;
		using Concs1DType = HDF5File::RaggedDataSet2D<ConcType>::Ragged2DType;

		/**
		 * Construct the group name for the given time step.
		 *
		 * @param concGroup The parent concentration group.
		 * @param loop The loop number in which it is.
		 * @param timeStep The time step the group will represent.
		 * @return A string to use for the name of the time step group for
		 *          the given time step.
		 */
		static std::string
		makeGroupName(const ConcentrationGroup& concGroup, int ctrlStep,
			int loop, int timeStep);

		/**
		 * Construct a TimestepGroup.
		 * Default and copy constructors explicitly disallowed.
		 */
		TimestepGroup(void) = delete;
		TimestepGroup(const TimestepGroup& other) = delete;

		/**
		 * Create and populate a Timestep group within the given
		 * concentration group.
		 *
		 * @param concGroup The concentration group
		 * @param loop The loop number
		 * @param timeStep The number of the time step
		 * @param time The physical time at this time step
		 * @param previousTime The physical time at the previous time step
		 * @param deltaTime The physical length of the time step
		 */
		TimestepGroup(const ConcentrationGroup& concGroup, int ctrlStep,
			int loop, int timeStep, double time, double previousTime,
			double deltaTime);

		/**
		 * Open a TimestepGroup within the given concentration group.
		 *
		 * @param concGroup The concentration group
		 * @param loop The loop number
		 * @param timeStep The time step of the desired group.
		 */
		TimestepGroup(const ConcentrationGroup& concGroup, int ctrlStep,
			int loop, int timeStep);

		/**
		 * Update a Timestep group within the given
		 * concentration group.
		 *
		 * @param time The physical time at this time step
		 * @param previousTime The physical time at the previous time step
		 * @param deltaTime The physical length of the time step
		 */
		void
		updateTimestepGroup(double time, double previousTime, double deltaTime);

		/**
		 * Save the grid information to our timestep group.
		 *
		 * @param grid The grid points in the x direction (depth)
		 * @param ny The number of grid points in the y direction
		 * @param hy The step size in the y direction
		 * @param nz The number of grid points in the z direction
		 * @param hz The step size in the z direction
		 */
		void
		writeGrid(const std::vector<double>& grid, int ny = 0, double hy = 0.0,
			int nz = 0, double hz = 0.0) const;

		/**
		 * Save the fluence information to our timestep group.
		 *
		 * @param fluence The vector of fluences
		 */
		void
		writeFluence(const std::vector<double>& fluence) const;

		/**
		 * Save the surface positions to our timestep group.
		 *
		 * @param nAtoms The quantity of atoms at the surface
		 * @param previousFluxes The previous fluxes
		 * @param atomNames The names for the atom types
		 */
		void
		writeSurface1D(std::vector<Data1DType> nAtoms,
			std::vector<Data1DType> previousFluxes,
			std::vector<std::string> atomNames) const;

		/**
		 * Save the surface positions to our timestep group.
		 *
		 * @param iSurface The indices of the surface position
		 * @param nAtoms The quantity of atoms at the surface
		 * @param previousFluxes The previous fluxes
		 * @param atomNames The names for the atom types
		 */
		void
		writeSurface2D(const Surface2DType& iSurface,
			std::vector<Data2DType> nAtoms,
			std::vector<Data2DType> previousFluxes,
			std::vector<std::string> atomNames) const;

		/**
		 * Save the surface positions to our timestep group.
		 *
		 * @param iSurface The indices of the surface position
		 * @param nAtoms The quantity of atoms at the surface
		 * @param previousFluxes The previous fluxes
		 * @param atomNames The names for the atom types
		 */
		void
		writeSurface3D(const Surface3DType& iSurface,
			std::vector<Data3DType> nAtoms,
			std::vector<Data3DType> previousFluxes,
			std::vector<std::string> atomNames) const;

		/**
		 * Save the bottom informations to our timestep group.
		 *
		 * @param nAtoms The quantity of atoms at the bulk
		 * @param previousFluxes The previous fluxes
		 * @param atomNames The names for the atom types
		 */
		void
		writeBottom1D(std::vector<Data1DType> nAtoms,
			std::vector<Data1DType> previousFluxes,
			std::vector<std::string> atomNames);

		/**
		 * Save the bottom informations to our timestep group.
		 *
		 * @param nAtoms The quantity of atoms at the bulk
		 * @param previousFluxes The previous fluxes
		 * @param atomNames The names for the atom types
		 */
		void
		writeBottom2D(std::vector<Data2DType> nAtoms,
			std::vector<Data2DType> previousFluxes,
			std::vector<std::string> atomNames);

		/**
		 * Save the bottom informations to our timestep group.
		 *
		 * @param nAtoms The quantity of atoms at the bulk
		 * @param previousFluxes The previous fluxes
		 * @param atomNames The names for the atom types
		 */
		void
		writeBottom3D(std::vector<Data3DType> nAtoms,
			std::vector<Data3DType> previousFluxes,
			std::vector<std::string> atomNames);

		/**
		 * Save the bursting informations to our timestep group.
		 *
		 * @param nHe The quantity of helium lost from bursting
		 * @param nD The quantity of deuterium lost from bursting
		 * @param nT The quantity of tritium lost from bursting
		 */
		void
		writeBursting(Data1DType nHe, Data1DType nD, Data1DType nT);

		/**
		 * Add a concentration dataset at a specific grid point.
		 *
		 * @param size The size of the dataset to create
		 * @param concArray The array of concentration at a grid point
		 * @param write To know if we own the data to write
		 * @param i The index of the position on the grid on the x direction
		 * @param j The index of the position on the grid on the y direction
		 * @param k The index of the position on the grid on the z direction
		 */
		// TODO this should go away.
		void
		writeConcentrationDataset(int size, double concArray[][2], bool write,
			int i, int j = -1, int k = -1);

		/**
		 * Add a concentration dataset for all grid points in a 1D problem.
		 * Caller gives us a 2D ragged representation, and we flatten
		 * it into a 1D dataset and add a 1D "starting index" array.
		 * Assumes that grid point slabs are assigned to processes in
		 * MPI rank order.
		 *
		 * @param file The HDF5 file that owns our group.  Needed to support
		 *              parallel file access.
		 * @param baseX Index of first grid point we own.
		 * @param concs Concentrations associated with grid points we own.
		 *              Must have size equal to number of grid points we own.
		 *              Element i contains concentration data for
		 *              (baseX + i)
		 */
		// TODO measure performance gain when caller gives us
		// flattened array instead of having us convert to/from flat
		// representation.
		void
		writeConcentrations(
			const XFile& file, int baseX, const Concs1DType& concs) const;

		/**
		 * Read concentration dataset for our grid points in a 1D problem.
		 * Assumes that grid point slabs are assigned to processes in
		 * MPI rank order.
		 *
		 * @param file The HDF5 file that owns our group.  Needed to support
		 *              parallel file access.
		 * @param baseX Index of first grid point we own.
		 * @param numX Number of grid points we own.
		 * @return Concentrations associated with grid points we own.
		 *              Element i contains concentration data for
		 *              (baseX + i)
		 */
		Concs1DType
		readConcentrations(const XFile& file, int baseX, int numX) const;

		/**
		 * Read the times from our timestep group.
		 *
		 * @return pair(time, deltaTime) containing the physical time to
		 *          be changed and the time step length to be changed.
		 */
		std::pair<double, double>
		readTimes(void) const;

		/**
		 * Read the previous time from our concentration group.
		 *
		 * @return The physical time at the previous timestep
		 */
		double
		readPreviousTime(void) const;

		/**
		 * Read the grid size.
		 *
		 * @param nx The number of grid points in the x direction (depth)
		 * @param hx The step size in the x direction
		 * @param ny The number of grid points in the y direction
		 * @param hy The step size in the y direction
		 * @param nz The number of grid points in the z direction
		 * @param hz The step size in the z direction
		 */
		void
		readSizes(int& nx, double& hx, int& ny, double& hy, int& nz,
			double& hz) const;

		/**
		 * Read the grid.
		 *
		 * @return The grid
		 */
		std::vector<double>
		readGrid() const;

		/**
		 * Read the fluences.
		 *
		 * @return The fluences
		 */
		std::vector<double>
		readFluence() const;

		/**
		 * Read the surface position from our concentration group in
		 * the case of a 2D grid (a vector of surface positions).
		 *
		 * @return The vector of indices of the surface position
		 */
		Surface2DType
		readSurface2D(void) const;

		/**
		 * Read the surface position from our concentration group in
		 * the case of a 3D grid (a vector of vector of surface positions).
		 *
		 * @return The vector of vector of indices of the surface position
		 */
		Surface3DType
		readSurface3D(void) const;

		/**
		 * Read some data from our concentration
		 * group in the case of a 1D grid (one float).
		 *
		 * @param dataName The name of the data we want
		 * @return The value of the data
		 */
		Data1DType
		readData1D(const std::string& dataName) const;

		/**
		 * Read some data from our concentration group in
		 * the case of a 2D grid (a vector).
		 *
		 * @param dataName The name of the data we want
		 * @return The vector of the data
		 */
		Data2DType
		readData2D(const std::string& dataName) const;

		/**
		 * Read some data from our concentration group file in
		 * the case of a 3D grid (a vector of vector).
		 *
		 * @param dataName The name of the data we want
		 * @return The vector of vector of data
		 */
		Data3DType
		readData3D(const std::string& dataName) const;

		/**
		 * Read our (i,j,k)-th grid point concentrations.
		 *
		 * @param i The index of the grid point on the x axis
		 * @param j The index of the grid point on the y axis
		 * @param k The index of the grid point on the z axis
		 * @return The vector of concentrations
		 */
		// TODO remove once have added support for 0D, 2D, and 3D
		// parallel reads of concentrations.
		Data3DType
		readGridPoint(int i, int j = -1, int k = -1) const;
	};

	// Our concentrations group.
	class ConcentrationGroup : public HDF5File::Group
	{
	private:
		// Name of our last timestep and loop attribute.
		static const std::string lastTimestepAttrName;
		static const std::string lastLoopAttrName;
		static const std::string lastCtrlStepAttrName;

	public:
		// Path of the concentrations group within the file.
		static const fs::path path;

		// Create or open the concentrationsGroup.
		ConcentrationGroup(void) = delete;
		ConcentrationGroup(const ConcentrationGroup& other) = delete;
		ConcentrationGroup(const XFile& file, bool create = false);

		/**
		 * Add a concentration timestep group for the given time step.
		 *
		 * @param loop The loop number
		 * @param timeStep The number of the time step
		 * @param time The physical time at this time step
		 * @param previousTime The physical time at the previous time step
		 * @param deltaTime The physical length of the time step
		 */
		std::unique_ptr<TimestepGroup>
		addTimestepGroup(int ctrlStep, int loop, int timeStep, double time,
			double previousTime, double deltaTime) const;

		/**
		 * Obtain the last control step known to our group.
		 *
		 * @return Control step of last TimestepGroup written to our group.
		 */
		int
		getLastControlStep(void) const;

		/**
		 * Obtain the last loop known to our group.
		 *
		 * @return Loop of last TimestepGroup written to our group.
		 */
		int
		getLastLoop(void) const;

		/**
		 * Obtain the last timestep known to our group.
		 *
		 * @return Time step of last TimestepGroup written to our group.
		 */
		int
		getLastTimeStep(void) const;

		/**
		 * Determine if we have any TimestepGroups.
		 *
		 * @return True if any TimestepGroups have been written.
		 */
		bool
		hasTimesteps(void) const
		{
			return getLastTimeStep() >= 0;
		}

		/**
		 * Access the TimestepGroup associated with the given time step.
		 *
		 * @param loop The loop number
		 * @param timeStep Time step of the desired TimestepGroup.
		 * @return TimestepGroup associated with the given time step.  Empty
		 *          pointer if the given time step is not known to us.
		 */
		std::unique_ptr<TimestepGroup>
		getTimestepGroup(int loop, int timeStep) const;

		/**
		 * Access the TimestepGroup associated with the last known time step.
		 *
		 * @return TimestepGroup associated with the last known time step.
		 *          Empty pointer if we do not yet have any time steps.
		 */
		std::unique_ptr<TimestepGroup>
		getLastTimestepGroup(void) const;
	};

	// A group describing a network within our HDF5 file.
	class NetworkGroup : public HDF5File::Group
	{
	public:
		// Concise name for type of network bounds
		// in HDF5 class method parameters.
		using NetworkBoundsType = std::vector<
			std::vector<core::network::IReactionNetwork::AmountType>>;

	private:
		// Names of network attribute.
		static const std::string sizeAttrName;
		static const std::string phaseSpaceAttrName;

	public:
		// Path to the network group within our HDF5 file.
		static const fs::path path;

		NetworkGroup(void) = delete;
		NetworkGroup(const NetworkGroup& other) = delete;

		/**
		 * Open an existing network group.
		 *
		 * @param file The file whose network group to open.
		 */
		NetworkGroup(const XFile& file);

		/**
		 * Creating a new network group.
		 *
		 * @param file The file where to create the network group.
		 * @param network The network to write.
		 */
		NetworkGroup(
			const XFile& file, core::network::IReactionNetwork& network);

		/**
		 * Read the network sizes from our group.
		 *
		 * @return The total size
		 */
		int
		readNetworkSize() const;

		/**
		 * Read the reactions for every cluster.
		 *
		 * @param network The network that need the reactions.
		 */
		void
		readReactions(core::network::IReactionNetwork& network) const;

		/**
		 * Copy ourself to the given file.
		 * A NetworkGroup must not already exist in the file.
		 *
		 * @param target The file to copy ourself to.
		 */
		void
		copyTo(const XFile& target) const;
	};

	// A group describing a cluster within our HDF5 file.
	class ClusterGroup : public HDF5File::Group
	{
	public:
		// Concise name for type of network bounds
		// in HDF5 class method parameters.
		using ClusterBoundsType =
			std::vector<core::network::IReactionNetwork::AmountType>;

	private:
		// Names of cluster attributes.
		static const std::string formationEnergyAttrName;
		static const std::string migrationEnergyAttrName;
		static const std::string diffusionFactorAttrName;
		static const std::string boundsAttrName;

	public:
		ClusterGroup(void) = delete;
		ClusterGroup(const ClusterGroup& other) = delete;

		/**
		 * Open a cluster group.
		 *
		 * @param networkGroup The group in which the cluster is located.
		 * @param id The id of the cluster.
		 */
		ClusterGroup(const NetworkGroup& networkGroup, int id);

		/**
		 * Create a cluster group.
		 *
		 * @param networkGroup The group in which to write.
		 * @param id The Cluster id
		 * @param bounds The bounds on the region covered by the cluster.
		 * @param formation The formation energy.
		 * @param migration The migration energy.
		 * @param diffusion The diffusion factor.
		 */
		ClusterGroup(const NetworkGroup& networkGroup, int id,
			ClusterBoundsType bounds, double formation, double migration,
			double diffusion);

		/**
		 * Construct the group name for the given time step.
		 *
		 * @param id The id of the cluster.
		 * @return A string to use for the name of the cluster group.
		 */
		static std::string
		makeGroupName(int id);

		/**
		 * Read the cluster properties from our group.
		 *
		 * @param formationEnergy The formation energy.
		 * @param migrationEnergy The migration energy.
		 * @param diffusionFactor The diffusion factor.
		 * @return The cluster bounds.
		 */
		ClusterBoundsType
		readCluster(double& formationEnergy, double& migrationEnergy,
			double& diffusionFactor) const;
	};

private:
	/**
	 * Pass through only Create* access modes.
	 *
	 * @param mode The access mode to check.
	 * @return The given access mode if it is a Create* mode.  Otherwise,
	 *          throws an exception.
	 */
	static AccessMode
	EnsureCreateAccessMode(AccessMode mode);

	/**
	 * Pass through only Open* access modes.
	 *
	 * @param mode The access mode to check.
	 * @return The given access mode if it is a Open* mode.  Otherwise,
	 *          throws an exception.
	 */
	static AccessMode
	EnsureOpenAccessMode(AccessMode mode);

public:
	/**
	 * Create and initialize a checkpoint file.
	 *
	 * @param path Path of file to create.
	 * @param create Dummy variable
	 * @param _comm The MPI communicator used to access the file.
	 * @param mode Access mode for file.  Only HDF5File Create* modes
	 *              are supported.
	 */
	XFile(fs::path path, int create, MPI_Comm _comm = MPI_COMM_WORLD,
		AccessMode mode = AccessMode::CreateOrTruncateIfExists);

	/**
	 * Open an existing checkpoint or network file.
	 *
	 * @param path Path of file to open.
	 * @param _comm The MPI communicator used to access the file.
	 * @param mode Access mode for file.  Only HDFFile Open* modes
	 *              are supported.
	 */
	XFile(fs::path path, MPI_Comm _comm = MPI_COMM_WORLD,
		AccessMode mode = AccessMode::OpenReadOnly);

	/**
	 * Check whether we have one of our top-level Groups.
	 *
	 * @return True iff we have the desired group.
	 */
	template <typename T>
	bool
	hasGroup(void) const
	{
		return HDF5File::hasGroup(T::path);
	}

	/**
	 * Access one of our top-level Groups within our file.
	 *
	 * @return The group object if we can open the group, else an empty pointer.
	 */
	template <typename T>
	std::unique_ptr<T>
	getGroup(void) const
	{
		std::unique_ptr<T> group;

		if (hasGroup<T>()) {
			// Open the group within our file.
			group = std::make_unique<T>(*this);
		}
		return std::move(group);
	}
};

} // namespace io
} // namespace xolotl

// Ensure we have definitions of template classes/methods.
#include <xolotl/io/XFileType.h>

#endif // XCORE_XFILE_H
