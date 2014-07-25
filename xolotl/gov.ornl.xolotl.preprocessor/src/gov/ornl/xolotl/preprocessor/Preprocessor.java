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
	private double[] heDiffusionFactors = { 0.0, 29.0, 32.0, 23.0, 17.0, 5.0,
			10.0 };

	// The migration energies for single species Helium clusters.
	private double[] heMigrationEnergies = { Double.POSITIVE_INFINITY, 0.13,
			0.20, 0.25, 0.20, 0.12, 0.3 };

	// The maximum size of a vacancy cluster in the network.
	private int maxV;

	// The diffusion factor for a single vacancy.
	private double vOneDiffusionFactor = 1800.0;

	// The migration energy for a single vacancy.
	private double vOneMigrationEnergy = 1.30;

	// The maximum size of an interstitial cluster in the network.
	private int maxI;

	// The maximum size of a mobile interstitial cluster.
	private int maxIDiffusionSize = 5;

	// The diffusion coefficients for single species interstitial clusters.
	private double[] iDiffusionFactors = { 0.0, 88.0, 80.0, 39.0, 20.0, 10.0 };

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
		Map<String, String> petscOptions = new HashMap<String, String>();
		petscOptions.put("-da_grid_x", "10");
		petscOptions.put("-ts_final_time", "1000");
		petscOptions.put("-ts_max_steps", "3");
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
		if ((petscList.get((petscList.size() - 1))).contains("-"))
			petscOptions.put(petscList.get((petscList.size() - 1)), "");

		// Loop through the Petsc list of strings to pair Petsc options with
		// their corresponding arguments and identify the stand-alone options
		for (int i = 1; i < petscList.size(); i++) {
			// Check if there is an option followed by a corresponding argument
			if (((petscList.get(i - 1)).contains("-"))
					&& (!(petscList.get(i)).contains("-"))) {
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
		if ( !(maxHe < 9) || !(maxHe > 0) ) {
			throw new IllegalArgumentException(
					"The maxium Helium size must be greater than 0 and less than 9 ( 0 < maxHeSize < 9 )");
		}

		// Set the maximum size of a vacancy cluster in the network.
		maxV = args.getMaxVSize();
		
		// The maximum size of an interstitial cluster in the network.
		maxI = args.getMaxISize();
		if ( !(maxI < 7) || !(maxI > 0) ) {
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
		if (args.isHeFluence())
			xolotlParams.setProperty("heFluence", args.getHeFluence());
		if (args.isCheckpoint())
			xolotlParams.setProperty("checkpoint", args.getCheckpoint());

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
	 * This operation generates the grid needed to write the concentrations.
	 * 
	 * @param dimension
	 *            The physical length of the grid
	 * @param refinement
	 *            The refinement of the grid
	 * 
	 * @return The array of physical positions on the grid
	 */
	public double[] generateGrid(int dimension, int refinement) {
		int totalLength = dimension;
		// Compute the total number of positions
		for (int i = 0; i < refinement; i++) {
			totalLength = (totalLength * 2) - 1;
		}
		// Create the array to return
		double[] toReturn = new double[totalLength];
		// Compute the distance between every position
		double increment = (double) dimension / (totalLength - 1);

		for (int i = 0; i < totalLength; i++) {
			toReturn[i] = i * increment;
		}

		return toReturn;
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
