/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.io.IOException;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Properties;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import uk.co.flamingpenguin.jewel.cli.*;
import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;

/**
 * This class generates a valid Xolotl input file with each line representing a
 * cluster with the following variables:
 * 
 * nHe nV nI E_He E_V E_I E_migration D_0
 * 
 * Each element is:
 * 
 * nHe - The number of Helium atoms in the cluster.
 * 
 * nV - The number of vacancies in the cluster.
 * 
 * nI - The number of interstitials in the cluster.
 * 
 * E_He - The binding energy of the cluster with He_1.
 * 
 * E_V - The binding energy of the cluster with V_1.
 * 
 * E_I - The binding energy of the cluster with I_1.
 * 
 * E_migration - The migration energy with which the cluster diffuses.
 * 
 * D_0 - The diffusion factor that sets the scale of the diffusion. (Used to
 * calculate the diffusion coefficient.)
 * 
 * All energies are in eV and all diffusion factors are in nm^2/s.
 * 
 * The preprocessor generates 2067 clusters with its default configuration. If
 * only the maximum number of vacancies is changed it generates a number of
 * clusters equal to
 * 
 * maxI + maxHe + 2053 + (maxV - 29) + 2(maxV - 29)(30 + maxV)
 * 
 * @author Jay Jay Billings
 * 
 */
public class Preprocessor {

	// The maximum size of a Helium cluster in the network.
	private int maxHe;

	// The maximum size of a mobile He cluster.
	private int maxHeDiffusionSize = 6;

	// The diffusion factors for single species Helium clusters.
	private double[] heDiffusionFactors = { 0.0, 2.9e+10, 3.2e+10, 2.3e+10,
			1.7e+10, 5.0e+09, 1.0e+09 };

	// The migration energies for single species Helium clusters.
	private double[] heMigrationEnergies = { Double.POSITIVE_INFINITY, 0.13,
			0.20, 0.25, 0.20, 0.12, 0.3 };

	// The maximum size of a vacancy cluster in the network.
	private int maxV;

	// The diffusion factor for a single vacancy.
	private double vOneDiffusionFactor = 1.8e+12;

	// The migration energy for a single vacancy.
	private double vOneMigrationEnergy = 1.30;

	// The maximum size of an interstitial cluster in the network.
	private int maxI;

	// The maximum size of a mobile interstitial cluster.
	private int maxIDiffusionSize = 5;

	// The diffusion coefficients for single species interstitial clusters.
	private double[] iDiffusionFactors = { 0.0, 8.8e+10, 8.0e+10, 3.9e+10,
			2.0e+10, 1.0e+10 };

	// The migration energies for single species interstitial clusters.
	private double[] iMigrationEnergies = { Double.POSITIVE_INFINITY, 0.01,
			0.02, 0.03, 0.04, 0.05 };

	/**
	 * The maximum number of Helium atoms that can be combined with a vacancy
	 * cluster with size equal to the index i in the array plus one. For
	 * example, an HeV size cluster with size 1 would have size = i+1 = 1 and i
	 * = 0. It could support a mixture of up to nine Helium atoms with one
	 * vacancy.
	 */
	private int[] maxHePerV = { 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60,
			65, 70, 75, 80, 85, 90, 95, 98, 100, 101, 103, 105, 107, 109, 110,
			112, 116 };

	/**
	 * The binding energy engine used to generate binding energies for the
	 * different clusters.
	 */
	private BindingEnergyEngine bindingEnergyEngine = new BindingEnergyEngine();

	/**
	 * The list of parameters that will be passed to Xolotl
	 */
	public Properties xolotlParams = new Properties();

	public Map<String, String> petscOptions = new HashMap<String, String>();

	/**
	 * This operation generates a single string containing the Petsc arguments
	 * that will be used to set the Petsc parameter required to run Xolotl.
	 * 
	 * @param petscArgs
	 *            The Petsc command line arguments as specified by the Arguments
	 *            interface
	 * @return A string containing the Petsc arguments that will be passed to
	 *         Xolotl
	 */
	private String generatePetscArgs(String petscArgs) {

		// Create a map of the default Petsc options and their corresponding
		// arguments, if any, where the key is the option and the value is
		// the argument
		petscOptions.put("-da_grid_x", "10");
		petscOptions.put("-ts_final_time", "50");
		petscOptions.put("-ts_dt", "1.0e-12");
		petscOptions.put("-ts_max_steps", "100");
		petscOptions.put("-ts_adapt_dt_max", "10");
		petscOptions.put("-ts_max_snes_failures", "200");
		petscOptions.put("-pc_type", "fieldsplit");
		petscOptions.put("-pc_fieldsplit_detect_coupling", "");
		petscOptions.put("-fieldsplit_0_pc_type", "redundant");
		petscOptions.put("-fieldsplit_1_pc_type", "sor");
		petscOptions.put("-snes_monitor", "");
		petscOptions.put("-ksp_monitor", "");
		petscOptions.put("-ts_monitor", "");

		// Get the string of Petsc arguments from the command line
		// and split the string around the blank spaces
		List<String> petscList = new ArrayList<String>();
		for (String str : petscArgs.split(" ")) {
			petscList.add(str);
		}

		// Check if the last string in the petscList is a stand-alone option
		String str = "-";
		char dash = str.charAt(0);
		if (((petscList.get((petscList.size() - 1))).charAt(0)) == dash) {
			petscOptions.put(petscList.get((petscList.size() - 1)), "");
		}

		// Loop through the Petsc list of strings to pair Petsc options with
		// their corresponding arguments and identify the stand-alone options
		for (int i = 1; i < petscList.size(); i++) {
			// Check if there is an option followed by a corresponding argument
			if (((petscList.get(i - 1)).charAt(0) == dash)
					&& !((petscList.get(i)).charAt(0) == dash)) {
				// Replace the default argument value with the specified value
				petscOptions.put(petscList.get(i - 1), petscList.get(i));
				i++;
			}
			// Identify stand-alone options
			else {
				petscOptions.put(petscList.get(i - 1), "");
			}
		}

		// Get the list of petscArgs and create a single string for them
		// in order to set the petsc parameter
		StringBuilder petscString = new StringBuilder();
		for (String key : petscOptions.keySet()) {
			if (petscOptions.get(key) == "")
				petscString.append(key + petscOptions.get(key) + " ");
			else
				petscString.append(key + " " + petscOptions.get(key) + " ");
		}

		return petscString.toString();
	}

	/**
	 * Constructor
	 * 
	 * @param args
	 *            The command line arguments as specified by the Arguments
	 *            interface.
	 */
	public Preprocessor(Arguments args) {

		// Set the maximum size of a Helium cluster in the network.
		maxHe = args.getMaxHeSize();
		// Check to make sure the user entered an appropriate value
		if (!(maxHe < 9) || !(maxHe > 0)) {
			throw new IllegalArgumentException(
					"The maxium Helium size must be greater than 0 and less than 9 ( 0 < maxHeSize < 9 )");
		}

		// Set the maximum size of a vacancy cluster in the network.
		maxV = args.getMaxVSize();

		// The maximum size of an interstitial cluster in the network.
		maxI = args.getMaxISize();
		if (!(maxI < 7) || !(maxI > 0)) {
			throw new IllegalArgumentException(
					"The maxium interstitial size must be greater than 0 and less than 7 ( 0 < maxHeSize < 7 )");
		}

		// Set the parameter options that will be passed to Xolotl
		xolotlParams.setProperty("startTemp", args.getStartTemp());
		xolotlParams.setProperty("networkFile", args.getNetworkFile());
		xolotlParams.setProperty("perfHandler", args.getPerfHandler());
		xolotlParams.setProperty("vizHandler", args.getVizHandler());
		xolotlParams.setProperty("petscArgs",
				generatePetscArgs(args.getPetscArgs()));

		// The following parameter options are optional and will only
		// be set if they are specified via the command line
		if (args.isMaterial())
			xolotlParams.setProperty("material", args.getMaterial());
		if (args.isTempFile())
			xolotlParams.setProperty("tempFile", args.getTempFile());
		if (args.isHeFlux())
			xolotlParams.setProperty("heFlux", args.getHeFlux());
		if (args.isMaxHeFluence())
			xolotlParams.setProperty("maxHeFluence", args.getMaxHeFluence());

	}

	/**
	 * This operation generates all Helium clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of Helium.
	 */
	private ArrayList<Cluster> generateHe() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Create the He clusters
		for (int i = 0; i < maxHe; i++) {
			// Create the cluster
			Cluster tmpCluster = new Cluster();
			tmpCluster.nHe = i + 1;
			tmpCluster.E_He = bindingEnergyEngine.getHeBindingEnergy(i + 1);
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Configure the diffusion parameters.
		for (int i = 0; i < maxHeDiffusionSize; i++) {
			Cluster tmpCluster = clusterList.get(i);
			tmpCluster.D_0 = heDiffusionFactors[i + 1];
			tmpCluster.E_m = heMigrationEnergies[i + 1];
		}

		return clusterList;
	}

	/**
	 * This operation creates an HeV cluster with the specified size. It
	 * configures the binding energies on its own.
	 * 
	 * @param heSize
	 *            The number of Helium atoms in the cluster
	 * @param vSize
	 *            The number of vacancies in the cluster
	 * @return The cluster.
	 */
	private Cluster makeHeVCluster(int heSize, int vSize) {
		// Create the cluster
		Cluster cluster = new Cluster();
		cluster.nHe = heSize;
		cluster.nV = vSize;
		// Treat everything like a mixed cluster and let the
		// BindingEnergyEngine delegate for single species clusters.
		cluster.E_He = bindingEnergyEngine.getHeVtoHeBindingEnergy(heSize,
				vSize);
		cluster.E_V = bindingEnergyEngine.getHeVtoVBindingEnergy(heSize, vSize);

		return cluster;
	}

	/**
	 * This operation generates all V and HeV clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of V and HeV clusters.
	 */
	private ArrayList<Cluster> generateHeV() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Loop over vacancies in the outer loop. Start at zero to account for
		// single He. This creates V and HeV up to the maximum size in the
		// maxHePerV array.
		for (int i = 1; i <= maxV && i <= maxHePerV.length; ++i) {
			// Add Helium
			for (int j = 0; j <= maxHePerV[i - 1]; j++) {
				// Add the cluster to the list
				clusterList.add(makeHeVCluster(j, i));
			}
		}

		// Create V and HeV up to the maximum length with a constant nHe/nV = 4.
		for (int i = maxHePerV.length + 1; i <= maxV; i++) {
			// Add Helium
			for (int j = 0; j <= i * 4; j++) {
				// Add the cluster to the list
				clusterList.add(makeHeVCluster(j, i));
			}
		}

		// Set V_1 diffusion parameters. V_1 is the first in the list, so it is
		// straightforward to set it.
		clusterList.get(0).D_0 = vOneDiffusionFactor;
		clusterList.get(0).E_m = vOneMigrationEnergy;

		return clusterList;
	};

	/**
	 * This operation generates all interstitials clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of interstitials.
	 */
	private ArrayList<Cluster> generateInterstitials() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Add the interstitials to the list
		for (int i = 0; i < maxI; ++i) {
			// Create the cluster
			Cluster tmpCluster = new Cluster();
			tmpCluster.nI = i + 1;
			tmpCluster.E_I = bindingEnergyEngine.getIBindingEnergy(i + 1);
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Configure the diffusion parameters.
		for (int i = 0; i < maxIDiffusionSize; i++) {
			Cluster tmpCluster = clusterList.get(i);
			tmpCluster.D_0 = iDiffusionFactors[i + 1];
			tmpCluster.E_m = iMigrationEnergies[i + 1];
		}

		return clusterList;
	};

	/**
	 * This operation creates an unstable cluster with the specified size.
	 * Unstable means that it will go through trap mutation with a relatively
	 * big dissociation rate.
	 * 
	 * @param heSize
	 *            The number of Helium atoms in the cluster
	 * @param vSize
	 *            The number of vacancies in the cluster
	 * @return The cluster.
	 */
	private Cluster makeUnstableCluster(int heSize, int vSize) {
		// Create the cluster
		Cluster cluster = new Cluster();
		cluster.nHe = heSize;
		cluster.nV = vSize;
		// Set its I binding energy to -5.0 to obtain a big dissociation rate
		// in Xolotl
		cluster.E_I = -5.0;

		return cluster;
	}

	/**
	 * This operation generates unstable clusters in the network. They are the
	 * ones just after the maximum number of He per V and will go under trap
	 * mutation. The last one should not be added to the network because it
	 * won't be able to dissociate into a cluster with a higher V number.
	 * 
	 * @return A list of clusters.
	 */
	private ArrayList<Cluster> generateUnstable() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Add the He one to the list
		clusterList.add(makeUnstableCluster(maxHe + 1, 0));

		// Loop over vacancies
		for (int i = 1; i < maxV && i <= maxHePerV.length; ++i) {
			// Add the unstable cluster to the list
			clusterList.add(makeUnstableCluster(maxHePerV[i - 1] + 1, i));
		}

		// Create unstable clusters for the cases where max He per V = 4.
		for (int i = maxHePerV.length + 1; i < maxV; i++) {
			// Add the unstable cluster to the list
			clusterList.add(makeUnstableCluster((i * 4) + 1, i));
		}

		return clusterList;
	};

	/**
	 * This operation generates the initial conditions based on the defaults and
	 * the incoming command line arguments.
	 * 
	 * @param args
	 *            The arguments read from the command line
	 * @return The list of clusters created by the preprocessor based on its
	 *         arguments and settings.
	 */
	public ArrayList<Cluster> generateNetwork(String[] args) {

		// Create the list of clusters
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Add all of the clusters to the list
		clusterList.addAll(generateInterstitials());
		clusterList.addAll(generateHe());
		clusterList.addAll(generateHeV());
		clusterList.addAll(generateUnstable());

		return clusterList;
	}

	/**
	 * This operation writes the parameters file that is needed to run Xolotl.
	 * 
	 * @param parameterFile
	 *            The parameter file name
	 * @param parameters
	 *            The parameters that will be written to the file
	 */
	public void writeParameterFile(String parameterFile, Properties parameters) {

		try {
			// Create the file containing the parameters
			FileOutputStream paramsFile = new FileOutputStream(parameterFile);

			// Write the parameters to the output file and save
			// the file to the project root folder
			parameters.store(paramsFile, null);
			// Flush the parameters to the intended stream
			paramsFile.flush();
			// Close the parameter file
			paramsFile.close();

		} catch (IOException io) {
			io.printStackTrace();
		}

		return;
	}

	/**
	 * This operation loads the parameters file that is needed to run Xolotl.
	 * 
	 * @param parameterFile
	 *            The parameter file name
	 */
	public Properties loadParameterFile(String parameterFile) {

		// Local declarations
		Properties inProperties = new Properties();

		try {

			FileInputStream inParamsFile = new FileInputStream(parameterFile);
			// Load the properties from the file
			inProperties.load(inParamsFile);
			// Close the parameter file
			inParamsFile.close();

		} catch (IOException io) {
			System.err.println("Error loading parameter file.");
			io.printStackTrace();
		}
		return inProperties;
	}

	/**
	 * This operation computes the total number of grid points
	 * 
	 * @param dimension
	 *            The physical length of the grid
	 * @param refinement
	 *            The refinement of the grid
	 * 
	 * @return The total number of points on the grid
	 */
	public int generateGrid(int dimension, int refinement) {
		int totalLength = dimension;
		// Compute the total number of positions
		for (int i = 0; i < refinement; i++) {
			totalLength = (totalLength * 2) - 1;
		}

		return totalLength;
	}

	/**
	 * This operation creates the HDF5 file needed by Xolotl
	 * 
	 * @param name
	 *            The name of the HDF5 file
	 */
	public void createHDF5(String name) {
		// The status of the previous HDF5 operation
		int status;

		try {
			// Create the HDF5 file
			int fileId = H5.H5Fcreate(name, HDF5Constants.H5F_ACC_TRUNC,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Close the HDF5 file
			status = H5.H5Fclose(fileId);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return;
	}

	/**
	 * This operation writes the header in the HDF5 file
	 * 
	 * @param name
	 *            The name of the HDF5 file
	 * @param dimension
	 *            The physical dimension of the grid
	 * @param refinement
	 *            The refinement of the grid
	 */
	public void writeHeader(String name, int[] dimension, int[] refinement) {
		// The status of the previous HDF5 operation
		int status;

		try {
			// Open the HDF5 file
			int fileId = H5.H5Fopen(name, HDF5Constants.H5F_ACC_RDWR,
					HDF5Constants.H5P_DEFAULT);

			// Create the header group
			int headerGroupId = H5.H5Gcreate(fileId, "headerGroup",
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);

			// Create, write, and close the physicalDim attribute
			int dimDataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int dimAttributeId = H5.H5Acreate(headerGroupId, "physicalDim",
					HDF5Constants.H5T_STD_I32LE, dimDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(dimAttributeId, HDF5Constants.H5T_STD_I32LE,
					dimension);
			status = H5.H5Aclose(dimAttributeId);

			// Create, write, and close the refinement attribute
			int refineDataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int refineAttributeId = H5.H5Acreate(headerGroupId, "refinement",
					HDF5Constants.H5T_STD_I32LE, refineDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(refineAttributeId,
					HDF5Constants.H5T_STD_I32LE, refinement);
			status = H5.H5Aclose(refineAttributeId);

			// Close everything
			status = H5.H5Gclose(headerGroupId);
			status = H5.H5Fclose(fileId);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return;
	}

	/**
	 * This operation copy the header from a given HDF5 file to the one created
	 * by the preprocessor
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which the header will be copied
	 * @param toName
	 *            The name of the created HDF5 file
	 * @return The total number of points of the grid
	 */
	public int copyHeader(String fromName, String toName) {
		// The status of the previous HDF5 operation
		int status;
		// The dimension and refinement
		int[] dimension = { 0 };
		int[] refinement = { 0 };

		try {
			// Read the header from the given file
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY,
					HDF5Constants.H5P_DEFAULT);

			// Open the header group
			int headerGroupId = H5.H5Gopen(fileId, "/headerGroup",
					HDF5Constants.H5P_DEFAULT);

			// Open and read the physicalDim attribute
			int dimAttributeId = H5.H5Aopen(headerGroupId, "physicalDim",
					HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(dimAttributeId, HDF5Constants.H5T_STD_I32LE,
					dimension);
			status = H5.H5Aclose(dimAttributeId);

			// Open and read the refinement attribute
			int refineAttributeId = H5.H5Aopen(headerGroupId, "refinement",
					HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(refineAttributeId, HDF5Constants.H5T_STD_I32LE,
					refinement);
			status = H5.H5Aclose(refineAttributeId);

			// Close everything
			status = H5.H5Gclose(headerGroupId);
			status = H5.H5Fclose(fileId);

			// Write the dimension and the refinement in the new HDF5 file
			// Open the HDF5 file
			fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR,
					HDF5Constants.H5P_DEFAULT);

			// Create the header group
			headerGroupId = H5.H5Gcreate(fileId, "headerGroup",
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);

			// Create, write, and close the physicalDim attribute
			int dimDataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			dimAttributeId = H5.H5Acreate(headerGroupId, "physicalDim",
					HDF5Constants.H5T_STD_I32LE, dimDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(dimAttributeId, HDF5Constants.H5T_STD_I32LE,
					dimension);
			status = H5.H5Aclose(dimAttributeId);

			// Create, write, and close the refinement attribute
			int refineDataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			refineAttributeId = H5.H5Acreate(headerGroupId, "refinement",
					HDF5Constants.H5T_STD_I32LE, refineDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(refineAttributeId,
					HDF5Constants.H5T_STD_I32LE, refinement);
			status = H5.H5Aclose(refineAttributeId);

			// Close everything
			status = H5.H5Gclose(headerGroupId);
			status = H5.H5Fclose(fileId);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return generateGrid(dimension[0], refinement[0]);
	}

	/**
	 * This operation read the value for the last time step saved in the given
	 * HDF5 file
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which concentrations
	 *            will be copied
	 * @return The value of the last time step saved in the HDF5 file
	 */
	public int readLastTimeStep(String fromName) {
		// The status of the previous HDF5 operation
		int status;
		// The indice of the last time step saved in the given HDF5 file
		int[] lastTimeStep = { -1 };

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY,
					HDF5Constants.H5P_DEFAULT);

			// Open the concentration group
			int concentrationGroupId = H5.H5Gopen(fileId,
					"/concentrationsGroup", HDF5Constants.H5P_DEFAULT);

			// Open and read the lastTimeStep attribute
			int lastAttributeId = H5.H5Aopen(concentrationGroupId,
					"lastTimeStep", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(lastAttributeId, HDF5Constants.H5T_STD_I32LE,
					lastTimeStep);
			status = H5.H5Aclose(lastAttributeId);

			// Close everything
			status = H5.H5Gclose(concentrationGroupId);
			status = H5.H5Fclose(fileId);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return lastTimeStep[0];
	}

	/**
	 * This operation creates the main concentration group in the new HDF5 file and 
	 * stores the last time step as being 0
	 * 
	 * @param toName
	 *            The name of the created HDF5 file
	 */
	public void createInitialConcentrationGroup(String toName) {
		// The status of the previous HDF5 operation
		int status;
		
		try {
			// Open the created HDF5 file to start writing the new
			// concentration group in it
			int fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR,
					HDF5Constants.H5P_DEFAULT);

			// Create the concentration group
			int newConcGroupId = H5.H5Gcreate(fileId,
					"concentrationsGroup", HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Create, write, and close the last written time step attribute
			// for the new file
			int[] newTimeStep = { 0 };
			int lastDataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int lastAttributeId = H5.H5Acreate(newConcGroupId, "lastTimeStep",
					HDF5Constants.H5T_STD_I32LE, lastDataspaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(lastAttributeId,
					HDF5Constants.H5T_STD_I32LE, newTimeStep);
			status = H5.H5Aclose(lastAttributeId);
			
			// Close everything
			status = H5.H5Gclose(newConcGroupId);
			status = H5.H5Fclose(fileId);
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return;
	}
	
	/**
	 * This operation reads the time attributes from the given HDF5 file
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which concentrations
	 *            will be copied
	 * @param lastTimeStep
	 *            The value of the last time step to know which 
	 *            concentration group to look at
	 * @return An array containing the absolute time and the timestep
	 */
	public double[] readTimes(String fromName, int lastTimeStep) {
		// The status of the previous HDF5 operation
		int status;
		// The array for the times
		double[] times = {0.0, 0.0};
		
		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY,
					HDF5Constants.H5P_DEFAULT);
			
			// Open the corresponding sub group
			// Set the name of the sub group
			String subGroupName = "concentrationsGroup/concentration_"
					+ lastTimeStep;

			// Open this specific concentration sub group
			int concentrationGroupId = H5.H5Gopen(fileId, subGroupName,
					HDF5Constants.H5P_DEFAULT);
			
			// Initialize the times
			double[] time = { 0.0 };
			double[] deltaTime = { 0.0 };

			// Open and read the absoluteTime attribute
			int timeAttributeId = H5.H5Aopen(concentrationGroupId,
					"absoluteTime", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(timeAttributeId,
					HDF5Constants.H5T_IEEE_F64LE, time);
			status = H5.H5Aclose(timeAttributeId);

			// Open and read the deltaTime attribute
			int deltaAttributeId = H5.H5Aopen(concentrationGroupId,
					"deltaTime", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(deltaAttributeId,
					HDF5Constants.H5T_IEEE_F64LE, deltaTime);
			status = H5.H5Aclose(deltaAttributeId);

			// Close everything
			status = H5.H5Gclose(concentrationGroupId);
			status = H5.H5Fclose(fileId);
			
			// Set the times in the array
			times[0] = time[0];
			times[1] = deltaTime[0];
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return times;
	}

	/**
	 * This operation creates the concentration subgroup in the new HDF5 file and 
	 * stores the absolute time and the time step
	 * 
	 * @param toName
	 *            The name of the created HDF5 file
	 * @param times
	 * 		      An array containing the absolute time and the timestep
	 */
	public void createConcentrationSubGroup(String toName, double[] times) {
		// The status of the previous HDF5 operation
		int status;
		
		try {
			// Open the created HDF5 file 
			int fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR,
					HDF5Constants.H5P_DEFAULT);
			
			// Create the concentration subgroup where the concentrations
			// will be actually stored
			int newConcSubGroupId = H5.H5Gcreate(fileId,
					"concentrationsGroup/concentration_0", HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			
			// Create, write, and close the absolute time attribute
			double[] time = { times[0] };
			int timeDataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int timeAttributeId = H5.H5Acreate(newConcSubGroupId,
					"absoluteTime", HDF5Constants.H5T_IEEE_F64LE,
					timeDataspaceId, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(timeAttributeId,
					HDF5Constants.H5T_IEEE_F64LE, time);
			status = H5.H5Sclose(timeDataspaceId);
			status = H5.H5Aclose(timeAttributeId);

			// Create, write, and close the timestep time attribute
			double[] deltaTime = { times[1] };
			int deltaDataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int deltaAttributeId = H5.H5Acreate(newConcSubGroupId, "deltaTime",
					HDF5Constants.H5T_IEEE_F64LE, deltaDataspaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(deltaAttributeId,
					HDF5Constants.H5T_IEEE_F64LE, deltaTime);
			status = H5.H5Sclose(deltaDataspaceId);
			status = H5.H5Aclose(deltaAttributeId);

			// Close everything
			status = H5.H5Gclose(newConcSubGroupId);
			status = H5.H5Fclose(fileId);
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return;
	}
	
	/**
	 * This operation reads the network of an already existing file and create 
	 * a map of the index and composition of the clusters 
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which the network will be read
	 * @return A map of the index and composition
	 */
	public int[][] readNetwork(String fromName) {
		// The status of the previous HDF5 operation
		int status;
		// The array for the times
		int[][] map = new int[0][0];
		
		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY,
					HDF5Constants.H5P_DEFAULT);

			// Open the dataset
			int datasetId = H5.H5Dopen(fileId, "/networkGroup/network", 
					HDF5Constants.H5P_DEFAULT);

			// Open and read the networkSize attribute
			int networkSizeAttributeId = H5.H5Aopen(datasetId, "networkSize", 
					HDF5Constants.H5P_DEFAULT);
			int[] networkSize = { 0 };
			status = H5.H5Aread(networkSizeAttributeId, HDF5Constants.H5T_STD_I32LE, networkSize);
			status = H5.H5Aclose(networkSizeAttributeId);

			// Create the array that will receive the network
			double[][] networkArray = new double[networkSize[0]][8];

			// Read the data set
			status = H5.H5Dread(datasetId, HDF5Constants.H5T_IEEE_F64LE, HDF5Constants.H5S_ALL, 
					HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, networkArray);

			// Close everything
			status = H5.H5Dclose(datasetId);
			status = H5.H5Fclose(fileId);
			
			// Recreate the map at the size of the network array
			map = new int[networkSize[0]][3];
			
			// Loop on the network array to fill the map
			for (int i = 0; i < networkSize[0]; i++) {
				map[i][0] = (int) networkArray[i][0]; // He
				map[i][1] = (int) networkArray[i][1]; // V
				map[i][2] = (int) networkArray[i][2]; // I
			}
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return map;
	}
	
	/**
	 * This operation checks whether a dataset exists for this position 
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which concentrations
	 *            will be copied
	 * @param lastTimeStep
	 *            The value of the last time step to know which 
	 *            concentration group to look at
	 * @param position
	 * 			  The index for the position
	 * @return True if there is a dataset at this position
	 */
	public boolean hasConcentrationDataset(String fromName, int lastTimeStep, int position) {
		// The status of the previous HDF5 operation
		int status;
		// Initializing the boolean
		boolean hasData = false;
		
		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY,
					HDF5Constants.H5P_DEFAULT);
			
			// Create the name of the dataset
			String datasetName = "concentrationsGroup/concentration_"
					+ lastTimeStep + "/position_" + position;
			// Check if the dataset exists
			hasData = H5.H5Lexists(fileId, datasetName,
					HDF5Constants.H5P_DEFAULT);

			// Close everything
			status = H5.H5Fclose(fileId);
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return hasData;
	}
	
	/**
	 * This operation reads the concentration at a specific position from 
	 * the given HDF5 file
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which concentrations
	 *            will be copied
	 * @param lastTimeStep
	 *            The value of the last time step to know which 
	 *            concentration group to look at
	 * @param position
	 * 			  The index for the position
	 * @return A two dimensional array containing the concentrations
	 */
	public double[][] readConcentration(String fromName, int lastTimeStep, int position) {
		// The status of the previous HDF5 operation
		int status;
		// The array for the times
		double[][] concentration = new double[0][0];
		
		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY,
					HDF5Constants.H5P_DEFAULT);
			
			// Create the name of the dataset
			String datasetName = "concentrationsGroup/concentration_"
					+ lastTimeStep + "/position_" + position;

			// Open it
			int datasetId = H5.H5Dopen(fileId, datasetName,
					HDF5Constants.H5P_DEFAULT);

			// Read the dataset length attribute
			int[] length = { -1 };
			int lengthAttributeId = H5.H5Aopen(datasetId,
					"datasetLength", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(lengthAttributeId,
					HDF5Constants.H5T_STD_I32LE, length);
			status = H5.H5Aclose(lengthAttributeId);

			// Create the array that will receive the concentrations
			concentration = new double[length[0]][2];

			// Read the data set
			status = H5.H5Dread(datasetId,
					HDF5Constants.H5T_IEEE_F64LE,
					HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL,
					HDF5Constants.H5P_DEFAULT, concentration);

			// Close everything
			status = H5.H5Dclose(datasetId);
			status = H5.H5Fclose(fileId);
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return concentration;
	}
	
	/**
	 * This operation writes the concentrations at this position in the new HDF5 file
	 * 
	 * @param toName
	 *            The name of the created HDF5 file
	 * @param position
	 * 			  The index for the position
	 * @param concentration
	 *            The two dimensional array containing the concentrations
	 * @param map
	 *            The map containing the index and composition of the clusters
	 *            in the previous network
	 * @param clusters
	 * 			  The list of clusters composing the new network
	 */
	public void writeConcentration(String toName, int position, 
			double[][] concentration, int[][] map, ArrayList<Cluster> clusters) {
		// The status of the previous HDF5 operation
		int status;
		// The list storing the new index and concentrations
		ArrayList<double[]> newConc = new ArrayList<double[]>();
		
		try {
			// Looping over the array on concentration to check that the clusters 
			// from the previous network still exist in the new network
			for (int i = 0; i < concentration.length; i++) {
				// boolean to know if the cluster is in the new network
				boolean found = false;
				// index in the new network
				int j = 0;
				// Loop on the new network
				for (Cluster cluster : clusters) {
					// Check the composition
					if ((cluster.nHe == map[(int) concentration[i][0]][0])
							&& (cluster.nV == map[(int) concentration[i][0]][1])
							&& (cluster.nI == map[(int) concentration[i][0]][2])) {
						// Add the cluster to the new list, with the new index 
						// but previous concentration
						double[] concCluster = { (double) j,  concentration[i][1] };
						newConc.add(concCluster);
						found = true;
					}
					// increment the j index
					j++;
				}
				
				// If the cluster was not found
				if (!found) {
					// Inform the user
					System.out.println("Cluster with the following composition (He, V, I): "
							+ map[(int) concentration[i][0]][0] + ", " 
							+ map[(int) concentration[i][0]][1] + ", " 
							+ map[(int) concentration[i][0]][2] 
							+ " is not present in the new network. Its concentration was "
							+ concentration[i][1] + " at the "
							+ position + "th grid point.");
				}
			}
			
			// Recreate the concentrattion array from the new concentration list
			concentration = new double[newConc.size()][2];
			// Loop on it to fill it
			for (int i = 0; i < concentration.length; i++) {
				double[] temp = newConc.get(i);
				concentration[i][0] = temp[0];
				concentration[i][1] = temp[1];
			}
			
			// Open the created HDF5 file 
			int fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR,
					HDF5Constants.H5P_DEFAULT);

			// Create the same dataset in the new file
			String datasetName = "concentrationsGroup/concentration_0/position_" + position;

			// Create the dataspace for the dataset with dimension dims
			long[] dims = new long[2];
			dims[0] = concentration.length;
			dims[1] = 2;
			int concDataspaceId = H5.H5Screate_simple(2, dims, null);

			// Create the dataset for the concentrations
			int datasetId = H5.H5Dcreate(fileId, datasetName,
					HDF5Constants.H5T_IEEE_F64LE, concDataspaceId,
					HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);

			// Write the concentration array in the dataset
			status = H5.H5Dwrite(datasetId,
					HDF5Constants.H5T_IEEE_F64LE,
					HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL,
					HDF5Constants.H5P_DEFAULT, concentration);

			// Create the attribute for the dataset length and write it
			int[] length = { concentration.length };
			int lengthDataSpaceId = H5
					.H5Screate(HDF5Constants.H5S_SCALAR);
			int lengthAttributeId = H5.H5Acreate(datasetId,
					"datasetLength", HDF5Constants.H5T_STD_I32LE,
					lengthDataSpaceId, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(lengthAttributeId,
					HDF5Constants.H5T_STD_I32LE, length);

			// Close everything
			status = H5.H5Sclose(lengthDataSpaceId);
			status = H5.H5Aclose(lengthAttributeId);
			status = H5.H5Sclose(concDataspaceId);
			status = H5.H5Dclose(datasetId);
			status = H5.H5Fclose(fileId);
			
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
		return;
	}

	/**
	 * This operation copy the concentrations from a given HDF5 file to the one
	 * created by the preprocessor
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which the concentrations
	 *             will be copied
	 * @param toName
	 *            The name of the created HDF5 file
	 * @param size
	 *            The total number of grid points
	 * @param clusters
	 *            The list of clusters that compose the network that will 
	 *            be written in the created HDF5 file
	 */
	public void copyConcentration(String fromName, String toName, int size, 
			ArrayList<Cluster> clusters) {
		// The status of the previous HDF5 operation
		int status;

		try {
			// Read the concentration group from the given file to get the value
			// of the last timestep
			int lastTimeStep = readLastTimeStep(fromName);

			if (lastTimeStep >= 0) {
				// Create the concentration group in the new HDF5 file
				createInitialConcentrationGroup(toName);

				double[] times = new double[2];
				
				// Read the times from the given HDF5 file
				times = readTimes(fromName, lastTimeStep);
				
				// Create the concentration subgroup in the new file and store
				// the times
				createConcentrationSubGroup(toName, times);
				
				// Create a map of the index and composition of the network present 
				// in the given HDF5 file
				int[][] networkMap = readNetwork(fromName);
				
				// Loop on all the position to read and copy the values of the
				// concentrations
				for (int i = 0; i < size; i++) {
					if (!hasConcentrationDataset(fromName, lastTimeStep, i))
						continue;
					
					// Read the concentrations at this position from the given HDF5 file
					double[][] conc = readConcentration(fromName, lastTimeStep, i);
					
					// Write the concentrations for this position in the new HDF5 file
					writeConcentration(toName, i, conc, networkMap, clusters);
				}
			}

			else {
				// Complain
				System.err
						.println("The given HDF5 file doesn't have stored concentration. Aborting.");
				return;
			}

		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return;
	}

	/**
	 * This operation writes the generated network in the HDF5 file
	 * 
	 * @param name
	 *            The name of the HDF5 file
	 * @param clusters
	 *            The list of clusters representing the network
	 */
	public void writeNetwork(String name, ArrayList<Cluster> clusters) {
		// The status of the previous HDF5 operation
		int status;

		try {
			// Open the HDF5 file
			int fileId = H5.H5Fopen(name, HDF5Constants.H5F_ACC_RDWR,
					HDF5Constants.H5P_DEFAULT);

			// Create the header group
			int networkGroupId = H5.H5Gcreate(fileId, "networkGroup",
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);

			// Create the array that will store the network
			int networkSize = clusters.size();
			double[][] networkArray = new double[networkSize][8];

			int id = 0;
			// Loop on the clusters
			for (Cluster cluster : clusters) {
				// Store the composition
				networkArray[id][0] = cluster.nHe;
				networkArray[id][1] = cluster.nV;
				networkArray[id][2] = cluster.nI;

				// Store the binding energies
				networkArray[id][3] = cluster.E_He;
				networkArray[id][4] = cluster.E_V;
				networkArray[id][5] = cluster.E_I;

				// Store the migration energy
				networkArray[id][6] = cluster.E_m;

				// Store the diffusion factor
				networkArray[id][7] = cluster.D_0;

				// increment the id number
				id++;
			}

			// Create the dataspace for the network with dimension dims
			long[] dims = new long[2];
			dims[0] = networkSize;
			dims[1] = 8;
			int networkDataSpaceId = H5.H5Screate_simple(2, dims, null);

			// Create the dataset for the network
			int datasetId = H5.H5Dcreate(networkGroupId, "network",
					HDF5Constants.H5T_IEEE_F64LE, networkDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);

			// Write networkArray in the dataset
			status = H5.H5Dwrite(datasetId, HDF5Constants.H5T_IEEE_F64LE,
					HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL,
					HDF5Constants.H5P_DEFAULT, networkArray);

			// Create the attribute for the network size
			int networkSizeDataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int networkSizeAttributeId = H5.H5Acreate(datasetId, "networkSize",
					HDF5Constants.H5T_STD_I32LE, networkSizeDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Write it
			int[] tempNetworkSize = { networkSize };
			status = H5.H5Awrite(networkSizeAttributeId,
					HDF5Constants.H5T_STD_I32LE, tempNetworkSize);

			// Close everything
			status = H5.H5Aclose(networkSizeAttributeId);
			status = H5.H5Dclose(datasetId);
			status = H5.H5Gclose(networkGroupId);
			status = H5.H5Fclose(fileId);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return;
	}
}
