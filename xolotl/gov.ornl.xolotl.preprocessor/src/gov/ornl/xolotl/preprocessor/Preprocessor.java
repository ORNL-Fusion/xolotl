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
 * nHe - The number of helium atoms in the cluster.
 * 
 * nXe - The number of xenon atoms in the cluster.
 * 
 * nV - The number of vacancies in the cluster.
 * 
 * nI - The number of interstitials in the cluster.
 * 
 * E_f - The formation energy of the cluster
 * 
 * E_migration - The migration energy with which the cluster diffuses.
 * 
 * D_0 - The diffusion factor that sets the scale of the diffusion. (Used to
 * calculate the diffusion coefficient.)
 * 
 * All energies are in eV and all diffusion factors are in nm^2/s.
 * 
 * The preprocessor generates 2067 clusters with its default configuration.
 * 
 * @author Jay Jay Billings
 * 
 */
public class Preprocessor {

	// The maximum size of a helium cluster in the network.
	private int maxHe;

	// The maximum size of a mobile He cluster.
	private int maxHeDiffusionSize = 7;

	// The diffusion factors for single species helium clusters.
	private double[] heDiffusionFactors = { 0.0, 2.9e+10, 3.2e+10, 2.3e+10, 1.7e+10, 5.0e+09, 1.0e+09, 5.0e+08 };

	// The migration energies for single species helium clusters.
	private double[] heMigrationEnergies = { Double.POSITIVE_INFINITY, 0.13, 0.20, 0.25, 0.20, 0.12, 0.3, 0.4 };

	// The maximum size of a xenon cluster in the network.
	private int maxXe;

	// The diffusion factor for a single xenon.
	private double xeOneDiffusionFactor = 5.0e-3;

	// The migration energy for a single xenon.
	private double xeOneMigrationEnergy = 0.0;

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
	private double[] iDiffusionFactors = { 0.0, 8.8e+10, 8.0e+10, 3.9e+10, 2.0e+10, 1.0e+10 };

	// The migration energies for single species interstitial clusters.
	private double[] iMigrationEnergies = { Double.POSITIVE_INFINITY, 0.01, 0.02, 0.03, 0.04, 0.05 };

	// Whether the phase-cut method is used for the network of not
	private boolean usePhaseCut = false;

	/**
	 * The maximum number of helium atoms that can be combined with a vacancy
	 * cluster with size equal to the index i in the array plus one. For
	 * example, an HeV size cluster with size 1 would have size = i+1 = 1 and i
	 * = 0. It could support a mixture of up to nine helium atoms with one
	 * vacancy.
	 */
	private int[] maxHePerV = { 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 98, 100, 101,
			103, 105, 107, 109, 110, 112, 116 };

	/**
	 * The formation energy engine used to generate formation energies for the
	 * different clusters.
	 */
	private FormationEnergyEngine formationEnergyEngine = new FormationEnergyEngine();

	/**
	 * The number of spacial dimensions.
	 */
	private int dim;

	/**
	 * The list of parameters that will be passed to Xolotl
	 */
	public Properties xolotlParams = new Properties();

	public Map<String, String> petscOptions = new HashMap<String, String>();

	/**
	 * This operation generates a single string containing the PETSc arguments
	 * that will be used to set the PETSc parameter required to run Xolotl.
	 * 
	 * @param petscArgs
	 *            The PETSc command line arguments as specified by the Arguments
	 *            interface
	 * @return A string containing the PETSc arguments that will be passed to
	 *         Xolotl
	 */
	private String generatePetscArgs(String petscArgs) {
		// Create a map of the default PETSc options and their corresponding
		// arguments, if any, where the key is the option and the value is
		// the argument
		petscOptions.put("-ts_final_time", "1.0");
		petscOptions.put("-ts_dt", "1.0e-12");
		petscOptions.put("-ts_max_steps", "100");
		petscOptions.put("-ts_adapt_dt_max", "1.0e-6");
		petscOptions.put("-ts_max_snes_failures", "200");
		petscOptions.put("-pc_type", "fieldsplit");
		petscOptions.put("-pc_fieldsplit_detect_coupling", "");
		petscOptions.put("-fieldsplit_0_pc_type", "sor");
		petscOptions.put("-fieldsplit_1_pc_type", "redundant");
		petscOptions.put("-ts_monitor", "");
		petscOptions.put("-ts_exact_final_time", "stepover");

		// Get the string of PETSc arguments from the command line
		// and split the string around the blank spaces
		List<String> petscList = new ArrayList<String>();
		for (String str : petscArgs.split(" ")) {
			petscList.add(str);
		}

		// Change the default preconditionner if we are not in 1D
		if (dim > 1) {
			petscList.add("-fieldsplit_1_pc_type");
			petscList.add("gamg");
			petscList.add("-fieldsplit_1_ksp_type");
			petscList.add("gmres");
			petscList.add("-ksp_type");
			petscList.add("fgmres");
			petscList.add("-fieldsplit_1_pc_gamg_sym_graph");
			petscList.add("true");
		}

		// Change the timesteps default for the NE problem
		if (maxXe > 0) {
			petscList.add("-ts_dt");
			petscList.add("1.0e-1");
			petscList.add("-ts_adapt_dt_max");
			petscList.add("5.0e5");
		}

		// Create the dash character
		String str = "-";
		char dash = str.charAt(0);

		// Loop through the PETSc list of strings to pair PETSc options with
		// their corresponding arguments and identify the stand-alone options
		for (int i = 0; i < petscList.size(); i++) {
			// Check that we are reading the name of the option
			// ("-" at the beginning)
			if (petscList.get(i).charAt(0) == dash) {
				// If there is another element after and the element doesn't
				// start with "-", this is an option with value
				if (i < petscList.size() - 1 && !((petscList.get(i + 1)).charAt(0) == dash)) {
					petscOptions.put(petscList.get(i), petscList.get(i + 1));

					// Skip the value now
					i++;
				}
				// Else the option is not associated to a value
				else {
					petscOptions.put(petscList.get(i), "");
				}
			}
		}

		// Get the list of petscArgs and create a single string for them
		// in order to set the petsc parameter
		StringBuilder petscString = new StringBuilder();
		for (String key : petscOptions.keySet()) {
			if (petscOptions.get(key).equals(""))
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
		// Set the maximum size of a helium cluster in the network.
		maxHe = args.getMaxHeSize();
		// Check to make sure the user entered an appropriate value
		if ((maxHe > 8) || (maxHe < 0)) {
			throw new IllegalArgumentException("The maxium helium size must be less than 9 ( 0 <= maxHeSize < 9 )");
		}

		// Set the maximum size of a vacancy cluster in the network.
		maxXe = args.getMaxXeSize();
		if (maxXe < 0) {
			throw new IllegalArgumentException("The maxium xenon size must be positive ( 0 <= maxXeSize )");
		}

		// Set the maximum size of a vacancy cluster in the network.
		maxV = args.getMaxVSize();
		if (maxV < 0) {
			throw new IllegalArgumentException("The maxium vacancy must be positive ( 0 <= maxVSize )");
		}

		// The maximum size of an interstitial cluster in the network.
		maxI = args.getMaxISize();
		if (maxI < 0) {
			throw new IllegalArgumentException("The maxium interstitial must be positive ( 0 <= maxISize )");
		}

		// Whether the phase-cut method will be used or not
		usePhaseCut = args.isPhaseCut();

		// The number of dimension for the problem to solve
		dim = Integer.parseInt(args.getDimensions());

		// Set the parameter options that will be passed to Xolotl
		xolotlParams.setProperty("dimensions", args.getDimensions());
		xolotlParams.setProperty("startTemp", args.getStartTemp());
		xolotlParams.setProperty("flux", args.getFlux());
		xolotlParams.setProperty("networkFile", args.getNetworkFile());
		xolotlParams.setProperty("material", args.getMaterial());
		xolotlParams.setProperty("perfHandler", args.getPerfHandler());
		xolotlParams.setProperty("vizHandler", args.getVizHandler());
		xolotlParams.setProperty("petscArgs", generatePetscArgs(args.getPetscArgs()));
		xolotlParams.setProperty("process", args.getProcess());

		// The following parameter options are optional and will only
		// be set if they are specified via the command line
		if (args.isTempFile())
			xolotlParams.setProperty("tempFile", args.getTempFile());
		if (args.isFluxFile())
			xolotlParams.setProperty("fluxFile", args.getTempFile());
		if (args.isInitialV())
			xolotlParams.setProperty("initialV", args.getInitialV());
		if (args.isVoidPortion())
			xolotlParams.setProperty("voidPortion", args.getVoidPortion());
		if (args.isRegularGrid())
			xolotlParams.setProperty("regularGrid", args.getRegularGrid());
		if (args.isGrain())
			xolotlParams.setProperty("grain", args.getGrain());
		if (args.isSputter())
			xolotlParams.setProperty("sputtering", args.getSputter());

	}

	/**
	 * This operation generates all helium clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of helium.
	 */
	private ArrayList<Cluster> generateHe() {
		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Create the He clusters
		for (int i = 0; i < maxHe; i++) {
			// Create the cluster
			Cluster tmpCluster = new Cluster();
			tmpCluster.nHe = i + 1;
			tmpCluster.E_f = formationEnergyEngine.getHeFormationEnergy(i + 1);
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Configure the diffusion parameters.
		for (int i = 0; i < Math.min(maxHeDiffusionSize, maxHe); i++) {
			Cluster tmpCluster = clusterList.get(i);
			tmpCluster.D_0 = heDiffusionFactors[i + 1];
			tmpCluster.E_m = heMigrationEnergies[i + 1];
		}

		return clusterList;
	}

	/**
	 * This operation generates all xenon clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of xenon.
	 */
	private ArrayList<Cluster> generateXe() {
		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Create the He clusters
		for (int i = 0; i < maxXe; i++) {
			// Create the cluster
			Cluster tmpCluster = new Cluster();
			tmpCluster.nXe = i + 1;
			tmpCluster.E_f = formationEnergyEngine.getXeFormationEnergy(i + 1);
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Set Xe_1 diffusion parameters. Xe_1 is the first in the list, so it
		// is
		// straightforward to set it.
		if (maxXe > 0) {
			clusterList.get(0).D_0 = xeOneDiffusionFactor;
			clusterList.get(0).E_m = xeOneMigrationEnergy;
		}

		return clusterList;
	}

	/**
	 * This operation creates an HeV cluster with the specified size. It
	 * configures the formation energies on its own.
	 * 
	 * @param heSize
	 *            The number of helium atoms in the cluster
	 * @param vSize
	 *            The number of vacancies in the cluster
	 * @return The cluster.
	 */
	private Cluster makeHeVCluster(int heSize, int vSize) {
		// Create the cluster
		Cluster cluster = new Cluster();
		cluster.nHe = heSize;
		cluster.nV = vSize;
		// Separate the case where it is simply a V cluster
		if (heSize == 0) {
			cluster.E_f = formationEnergyEngine.getVFormationEnergy(vSize);
		} else {
			cluster.E_f = formationEnergyEngine.getHeVFormationEnergy(heSize, vSize);
		}

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

		// Check if the phase-cut method need to be applied
		if (usePhaseCut) {
			// Loop over vacancies in the outer loop.
			// This creates V and HeV.
			for (int i = 1; i <= maxV && i <= maxHePerV.length; ++i) {
				// For V < 12 loop on all the possible helium numbers up to
				// the maximum size in the maxHePerV array
				if (i < 12) {
					for (int j = 0; j <= maxHePerV[i - 1]; j++) {
						// Add the cluster to the list
						clusterList.add(makeHeVCluster(j, i));
					}
				}
				// For bigger V only add the needed helium sizes
				else {
					// Add the vacancy cluster
					clusterList.add(makeHeVCluster(0, i));
					for (int j = maxHePerV[i - 2] + 1; j <= maxHePerV[i - 1]; j++) {
						// Add the cluster to the list
						clusterList.add(makeHeVCluster(j, i));
					}
				}
			}

			// Create V and HeV up to the maximum length with a constant nHe/nV
			// = 4,
			// Keeping only the last four ones
			for (int i = maxHePerV.length + 1; i <= maxV; i++) {
				// Add the vacancy cluster
				clusterList.add(makeHeVCluster(0, i));
				for (int j = (i * 4) - 3; j <= i * 4; j++) {
					// Add the cluster to the list
					clusterList.add(makeHeVCluster(j, i));
				}
			}
		}
		// Else use the full network
		else {
			// Loop over vacancies in the outer loop.
			// This creates V and HeV up to the maximum size in the
			// maxHePerV array.
			for (int i = 1; i <= maxV && i <= maxHePerV.length; ++i) {
				// Loop on the helium number
				for (int j = 0; j <= maxHePerV[i - 1]; j++) {
					// Add the cluster to the list
					clusterList.add(makeHeVCluster(j, i));
				}
			}

			// Create V and HeV up to the maximum length with a constant nHe/nV
			// = 4.
			for (int i = maxHePerV.length + 1; i <= maxV; i++) {
				// Loop on the helium number
				for (int j = 0; j <= i * 4; j++) {
					// Add the cluster to the list
					clusterList.add(makeHeVCluster(j, i));
				}
			}
		}

		// Set V_1 diffusion parameters. V_1 is the first in the list, so it is
		// straightforward to set it.
		if (maxV > 0) {
			clusterList.get(0).D_0 = vOneDiffusionFactor;
			clusterList.get(0).E_m = vOneMigrationEnergy;
		}

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
			tmpCluster.E_f = formationEnergyEngine.getIFormationEnergy(i + 1);
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Configure the diffusion parameters.
		for (int i = 0; i < Math.min(maxIDiffusionSize, maxI); i++) {
			Cluster tmpCluster = clusterList.get(i);
			tmpCluster.D_0 = iDiffusionFactors[i + 1];
			tmpCluster.E_m = iMigrationEnergies[i + 1];
		}

		return clusterList;
	}

	/**
	 * This operation generates the initial conditions based on the defaults and
	 * the incoming command line arguments.
	 * 
	 * @return The list of clusters created by the preprocessor based on its
	 *         arguments and settings.
	 */
	public ArrayList<Cluster> generateNetwork() {
		// Create the list of clusters
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Add all of the clusters to the list
		clusterList.addAll(generateInterstitials());
		clusterList.addAll(generateHe());
		clusterList.addAll(generateHeV());
		clusterList.addAll(generateXe());

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
	 * This operation creates the HDF5 file needed by Xolotl
	 * 
	 * @param name
	 *            The name of the HDF5 file
	 */
	public void createHDF5(String name) {
		try {
			// Create the HDF5 file
			int fileId = H5.H5Fcreate(name, HDF5Constants.H5F_ACC_TRUNC, HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT);

			// Close the HDF5 file
			int status = H5.H5Fclose(fileId);
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
	 * @param args
	 *            The command line arguments as specified by the Arguments
	 *            interface.
	 */
	public void writeHeader(String name, Arguments args) {
		// Get the grid size
		int[] nxGrid = { args.getNxGrid() };
		int[] nyGrid = { args.getNyGrid() };
		int[] nzGrid = { args.getNzGrid() };
		// Get the step sizes
		double[] hx = { args.getXStepSize() };
		double[] hy = { args.getYStepSize() };
		double[] hz = { args.getZStepSize() };

		try {
			// Open the HDF5 file
			int fileId = H5.H5Fopen(name, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);

			// Create the header group
			int headerGroupId = H5.H5Gcreate(fileId, "headerGroup", HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Create, write, and close the nx attribute (nxGrid)
			int dataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int attributeId = H5.H5Acreate(headerGroupId, "nx", HDF5Constants.H5T_STD_I32LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			int status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, nxGrid);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the hx attribute (h)
			attributeId = H5.H5Acreate(headerGroupId, "hx", HDF5Constants.H5T_IEEE_F64LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_IEEE_F64LE, hx);
			status = H5.H5Aclose(attributeId);

			// Create, write, and close the ny attribute (nyGrid)
			attributeId = H5.H5Acreate(headerGroupId, "ny", HDF5Constants.H5T_STD_I32LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, nyGrid);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the hy attribute (h)
			attributeId = H5.H5Acreate(headerGroupId, "hy", HDF5Constants.H5T_IEEE_F64LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_IEEE_F64LE, hy);
			status = H5.H5Aclose(attributeId);

			// Create, write, and close the nz attribute (nzGrid)
			attributeId = H5.H5Acreate(headerGroupId, "nz", HDF5Constants.H5T_STD_I32LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, nzGrid);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the hz attribute (h)
			attributeId = H5.H5Acreate(headerGroupId, "hz", HDF5Constants.H5T_IEEE_F64LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_IEEE_F64LE, hz);
			status = H5.H5Aclose(attributeId);

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
	 * @return The array of sizes of the grid
	 */
	public int[] copyHeader(String fromName, String toName) {
		// The number of grid points and the step size
		int[] n = { 0 };
		double[] h = { 0 };
		// The grid size to return
		int[] gridSize = { 0, 0, 0 };

		try {
			// Read the header from the given file
			// Open the file
			int fileFromId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Open the header group
			int headerFromGroupId = H5.H5Gopen(fileFromId, "/headerGroup", HDF5Constants.H5P_DEFAULT);

			// Open the file where we want to write the information
			int fileToId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);

			// Create the header group
			int headerToGroupId = H5.H5Gcreate(fileToId, "headerGroup", HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Open, read, and close the nx attribute
			int attributeId = H5.H5Aopen(headerFromGroupId, "nx", HDF5Constants.H5P_DEFAULT);
			int status = H5.H5Aread(attributeId, HDF5Constants.H5T_STD_I32LE, n);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the nx attribute
			int dataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			attributeId = H5.H5Acreate(headerToGroupId, "nx", HDF5Constants.H5T_STD_I32LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, n);
			status = H5.H5Aclose(attributeId);

			// Store the size of the grid in the x direction
			gridSize[0] = n[0];

			// Open, read, and close the hx attribute
			attributeId = H5.H5Aopen(headerFromGroupId, "hx", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(attributeId, HDF5Constants.H5T_IEEE_F64LE, h);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the nx attribute
			attributeId = H5.H5Acreate(headerToGroupId, "hx", HDF5Constants.H5T_IEEE_F64LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_IEEE_F64LE, h);
			status = H5.H5Aclose(attributeId);

			// Open, read, and close the ny attribute
			attributeId = H5.H5Aopen(headerFromGroupId, "ny", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(attributeId, HDF5Constants.H5T_STD_I32LE, n);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the ny attribute
			attributeId = H5.H5Acreate(headerToGroupId, "ny", HDF5Constants.H5T_STD_I32LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, n);
			status = H5.H5Aclose(attributeId);

			// Store the size of the grid in the y direction
			gridSize[1] = n[0];

			// Open, read, and close the hy attribute
			attributeId = H5.H5Aopen(headerFromGroupId, "hy", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(attributeId, HDF5Constants.H5T_IEEE_F64LE, h);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the ny attribute
			attributeId = H5.H5Acreate(headerToGroupId, "hy", HDF5Constants.H5T_IEEE_F64LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_IEEE_F64LE, h);
			status = H5.H5Aclose(attributeId);

			// Open, read, and close the nz attribute
			attributeId = H5.H5Aopen(headerFromGroupId, "nz", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(attributeId, HDF5Constants.H5T_STD_I32LE, n);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the nz attribute
			attributeId = H5.H5Acreate(headerToGroupId, "nz", HDF5Constants.H5T_STD_I32LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, n);
			status = H5.H5Aclose(attributeId);

			// Store the size of the grid in the y direction
			gridSize[2] = n[0];

			// Open, read, and close the hz attribute
			attributeId = H5.H5Aopen(headerFromGroupId, "hz", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(attributeId, HDF5Constants.H5T_IEEE_F64LE, h);
			status = H5.H5Aclose(attributeId);
			// Create, write, and close the nz attribute
			attributeId = H5.H5Acreate(headerToGroupId, "hz", HDF5Constants.H5T_IEEE_F64LE, dataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(attributeId, HDF5Constants.H5T_IEEE_F64LE, h);
			status = H5.H5Aclose(attributeId);

			// Close everything
			status = H5.H5Gclose(headerFromGroupId);
			status = H5.H5Fclose(fileFromId);
			status = H5.H5Gclose(headerToGroupId);
			status = H5.H5Fclose(fileToId);
		} catch (RuntimeException e) {
			// Complain
			e.printStackTrace();
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return gridSize;
	}

	/**
	 * This operation read the value for the last time step saved in the given
	 * HDF5 file
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which concentrations will be
	 *            copied
	 * @return The value of the last time step saved in the HDF5 file
	 */
	public int readLastTimeStep(String fromName) {
		// The indice of the last time step saved in the given HDF5 file
		int[] lastTimeStep = { -1 };

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Open the concentration group
			int concentrationGroupId = H5.H5Gopen(fileId, "/concentrationsGroup", HDF5Constants.H5P_DEFAULT);

			// Open and read the lastTimeStep attribute
			int lastAttributeId = H5.H5Aopen(concentrationGroupId, "lastTimeStep", HDF5Constants.H5P_DEFAULT);
			int status = H5.H5Aread(lastAttributeId, HDF5Constants.H5T_STD_I32LE, lastTimeStep);
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
	 * This operation creates the main concentration group in the new HDF5 file
	 * and stores the last time step as being 0
	 * 
	 * @param toName
	 *            The name of the created HDF5 file
	 */
	public void createInitialConcentrationGroup(String toName) {
		try {
			// Open the created HDF5 file to start writing the new
			// concentration group in it
			int fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);

			// Create the concentration group
			int newConcGroupId = H5.H5Gcreate(fileId, "concentrationsGroup", HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Create, write, and close the last written time step attribute
			// for the new file
			int[] newTimeStep = { 0 };
			int lastDataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int lastAttributeId = H5.H5Acreate(newConcGroupId, "lastTimeStep", HDF5Constants.H5T_STD_I32LE,
					lastDataspaceId, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			int status = H5.H5Awrite(lastAttributeId, HDF5Constants.H5T_STD_I32LE, newTimeStep);
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
	 *            The name of the HDF5 file from which concentrations will be
	 *            copied
	 * @param lastTimeStep
	 *            The value of the last time step to know which concentration
	 *            group to look at
	 * @return An array containing the absolute time and the timestep
	 */
	public double[] readTimes(String fromName, int lastTimeStep) {
		// The array for the times
		double[] times = { 0.0, 0.0 };

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Open the corresponding sub group
			// Set the name of the sub group
			String subGroupName = "concentrationsGroup/concentration_" + lastTimeStep;

			// Open this specific concentration sub group
			int concentrationGroupId = H5.H5Gopen(fileId, subGroupName, HDF5Constants.H5P_DEFAULT);

			// Initialize the times
			double[] time = { 0.0 };
			double[] deltaTime = { 0.0 };

			// Open and read the absoluteTime attribute
			int timeAttributeId = H5.H5Aopen(concentrationGroupId, "absoluteTime", HDF5Constants.H5P_DEFAULT);
			int status = H5.H5Aread(timeAttributeId, HDF5Constants.H5T_IEEE_F64LE, time);
			status = H5.H5Aclose(timeAttributeId);

			// Open and read the deltaTime attribute
			int deltaAttributeId = H5.H5Aopen(concentrationGroupId, "deltaTime", HDF5Constants.H5P_DEFAULT);
			status = H5.H5Aread(deltaAttributeId, HDF5Constants.H5T_IEEE_F64LE, deltaTime);
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
	 * This operation creates the concentration subgroup in the new HDF5 file
	 * and stores the absolute time and the time step
	 * 
	 * @param toName
	 *            The name of the created HDF5 file
	 * @param times
	 *            An array containing the absolute time and the timestep
	 */
	public void createConcentrationSubGroup(String toName, double[] times) {
		try {
			// Open the created HDF5 file
			int fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);

			// Create the concentration subgroup where the concentrations
			// will be actually stored
			int newConcSubGroupId = H5.H5Gcreate(fileId, "concentrationsGroup/concentration_0",
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Create, write, and close the absolute time attribute
			double[] time = { times[0] };
			int timeDataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int timeAttributeId = H5.H5Acreate(newConcSubGroupId, "absoluteTime", HDF5Constants.H5T_IEEE_F64LE,
					timeDataspaceId, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			int status = H5.H5Awrite(timeAttributeId, HDF5Constants.H5T_IEEE_F64LE, time);
			status = H5.H5Sclose(timeDataspaceId);
			status = H5.H5Aclose(timeAttributeId);

			// Create, write, and close the timestep time attribute
			double[] deltaTime = { times[1] };
			int deltaDataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int deltaAttributeId = H5.H5Acreate(newConcSubGroupId, "deltaTime", HDF5Constants.H5T_IEEE_F64LE,
					deltaDataspaceId, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
			status = H5.H5Awrite(deltaAttributeId, HDF5Constants.H5T_IEEE_F64LE, deltaTime);
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
	 * This operation copies the surface position information from the already
	 * existing file to the concentration subgroup of the newly created file
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which the surface position
	 *            information will be read
	 * @param lastTimeStep
	 *            The value of the last time step to know which concentration
	 *            group to look at
	 * @param toName
	 *            The name of the created HDF5 file
	 */
	public void copySurface(String fromName, int lastTimeStep, String toName) {
		// The status of the previous HDF5 operation
		int status;

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Open the corresponding sub group
			// Set the name of the sub group
			String subGroupName = "concentrationsGroup/concentration_" + lastTimeStep;

			// Open this specific concentration sub group
			int concentrationGroupId = H5.H5Gopen(fileId, subGroupName, HDF5Constants.H5P_DEFAULT);

			// Check if the surface position is an attribute
			// (thus the Xolotl was ran in 1D)
			if (H5.H5Aexists(concentrationGroupId, "iSurface")) {
				// Read the attribute
				int[] iSurface = { -1 };
				int attributeId = H5.H5Aopen(concentrationGroupId, "iSurface", HDF5Constants.H5P_DEFAULT);
				status = H5.H5Aread(attributeId, HDF5Constants.H5T_STD_I32LE, iSurface);

				// Close everything from the old file
				status = H5.H5Aclose(attributeId);
				status = H5.H5Gclose(concentrationGroupId);
				status = H5.H5Fclose(fileId);

				// And write the surface attribute in the new file
				// Open the created HDF5 file
				fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);
				// Open the new concentration subgroup
				concentrationGroupId = H5.H5Gopen(fileId, "concentrationsGroup/concentration_0",
						HDF5Constants.H5P_DEFAULT);
				// Create and write the surface position attribute
				int dataspaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
				attributeId = H5.H5Acreate(concentrationGroupId, "iSurface", HDF5Constants.H5T_STD_I32LE, dataspaceId,
						HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
				status = H5.H5Awrite(attributeId, HDF5Constants.H5T_STD_I32LE, iSurface);

				// Close everything from the new file
				status = H5.H5Aclose(attributeId);
				status = H5.H5Sclose(dataspaceId);
				status = H5.H5Gclose(concentrationGroupId);
				status = H5.H5Fclose(fileId);
			}
			// If the surface position information is not an attribute, it is a
			// dataset
			else {
				// Open the dataset
				int datasetId = H5.H5Dopen(concentrationGroupId, "iSurface", HDF5Constants.H5P_DEFAULT);
				// Get the dataspace from it
				int dataspaceId = H5.H5Dget_space(datasetId);
				// Determine the dimensionality
				int n = H5.H5Sget_simple_extent_ndims(dataspaceId);

				// Two possibilities now: either the dataset is a one dimension
				// one or
				// a two dimensions one, the former case corresponding to Xolotl
				// running in 2D,
				// the latter in 3D
				if (n == 1) {
					// We are in the Xolotl 2D case
					// Initialize the dimensions of the dataset
					long[] dims = { 0 };
					// Get the dimension of the dataset
					status = H5.H5Sget_simple_extent_dims(dataspaceId, dims, null);
					// Create the array that will receive the surface positions
					// information
					int[] surfaceArray = new int[(int) dims[0]];

					// Read the dataset
					status = H5.H5Dread(datasetId, HDF5Constants.H5T_STD_I32LE, HDF5Constants.H5S_ALL,
							HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, surfaceArray);

					// Close everything from the old file
					status = H5.H5Dclose(datasetId);
					status = H5.H5Sclose(dataspaceId);
					status = H5.H5Gclose(concentrationGroupId);
					status = H5.H5Fclose(fileId);

					// Write this surface position dataset into the new file
					// Open the created HDF5 file
					fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);
					// Open the new concentration subgroup
					concentrationGroupId = H5.H5Gopen(fileId, "concentrationsGroup/concentration_0",
							HDF5Constants.H5P_DEFAULT);
					// Create the dataspace
					dataspaceId = H5.H5Screate_simple(1, dims, null);
					// Create the dataset for the concentrations
					datasetId = H5.H5Dcreate(concentrationGroupId, "iSurface", HDF5Constants.H5T_STD_I32LE, dataspaceId,
							HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
					// Write the surface position array in the dataset
					status = H5.H5Dwrite(datasetId, HDF5Constants.H5T_STD_I32LE, HDF5Constants.H5S_ALL,
							HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, surfaceArray);

					// Close everything from the new file
					status = H5.H5Dclose(datasetId);
					status = H5.H5Sclose(dataspaceId);
					status = H5.H5Gclose(concentrationGroupId);
					status = H5.H5Fclose(fileId);
				} else if (n == 2) {
					// We are in the Xolotl 3D case
					// Initialize the dimensions of the dataset
					long[] dims = { 0, 0 };
					// Get the dimension of the dataset
					status = H5.H5Sget_simple_extent_dims(dataspaceId, dims, null);
					// Create the array that will receive the surface positions
					// information
					int[][] surfaceArray = new int[(int) dims[0]][(int) dims[1]];

					// Read the dataset
					status = H5.H5Dread(datasetId, HDF5Constants.H5T_STD_I32LE, HDF5Constants.H5S_ALL,
							HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, surfaceArray);

					// Close everything from the old file
					status = H5.H5Dclose(datasetId);
					status = H5.H5Sclose(dataspaceId);
					status = H5.H5Gclose(concentrationGroupId);
					status = H5.H5Fclose(fileId);

					// Write this surface position dataset into the new file
					// Open the created HDF5 file
					fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);
					// Open the new concentration subgroup
					concentrationGroupId = H5.H5Gopen(fileId, "concentrationsGroup/concentration_0",
							HDF5Constants.H5P_DEFAULT);
					// Create the dataspace
					dataspaceId = H5.H5Screate_simple(2, dims, null);
					// Create the dataset for the concentrations
					datasetId = H5.H5Dcreate(concentrationGroupId, "iSurface", HDF5Constants.H5T_STD_I32LE, dataspaceId,
							HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
					// Write the surface position array in the dataset
					status = H5.H5Dwrite(datasetId, HDF5Constants.H5T_STD_I32LE, HDF5Constants.H5S_ALL,
							HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, surfaceArray);

					// Close everything from the new file
					status = H5.H5Dclose(datasetId);
					status = H5.H5Sclose(dataspaceId);
					status = H5.H5Gclose(concentrationGroupId);
					status = H5.H5Fclose(fileId);
				}
			}

		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return;

	}

	/**
	 * This operation reads the network of an already existing file and create a
	 * map of the index and composition of the clusters
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which the network will be read
	 * @return A map of the index and composition
	 */
	public int[][] readNetwork(String fromName) {
		// The array for the times
		int[][] map = new int[0][0];

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Open the dataset
			int datasetId = H5.H5Dopen(fileId, "/networkGroup/network", HDF5Constants.H5P_DEFAULT);

			// Open and read the networkSize attribute
			int networkSizeAttributeId = H5.H5Aopen(datasetId, "networkSize", HDF5Constants.H5P_DEFAULT);
			int[] networkSize = { 0 };
			int status = H5.H5Aread(networkSizeAttributeId, HDF5Constants.H5T_STD_I32LE, networkSize);
			status = H5.H5Aclose(networkSizeAttributeId);

			// Create the array that will receive the network
			double[][] networkArray = new double[networkSize[0]][6];

			// Read the data set
			status = H5.H5Dread(datasetId, HDF5Constants.H5T_IEEE_F64LE, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL,
					HDF5Constants.H5P_DEFAULT, networkArray);

			// Close everything
			status = H5.H5Dclose(datasetId);
			status = H5.H5Fclose(fileId);

			// Recreate the map at the size of the network array
			map = new int[networkSize[0]][3];

			// Loop on the network array to fill the map
			for (int i = 0; i < networkSize[0]; i++) {
				map[i][0] = (int) networkArray[i][0]; // He or Xe
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
	 *            The name of the HDF5 file from which concentrations will be
	 *            copied
	 * @param lastTimeStep
	 *            The value of the last time step to know which concentration
	 *            group to look at
	 * @param i
	 *            The index for the position on the x axis
	 * @param j
	 *            The index for the position on the y axis
	 * @param k
	 *            The index for the position on the z axis
	 * @return True if there is a dataset at this position
	 */
	public boolean hasConcentrationDataset(String fromName, int lastTimeStep, int i, int j, int k) {
		// Initializing the boolean
		boolean hasData = false;

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Create the name of the dataset
			String datasetName = "concentrationsGroup/concentration_" + lastTimeStep + "/position_" + i + "_" + j + "_"
					+ k;
			// Check if the dataset exists
			hasData = H5.H5Lexists(fileId, datasetName, HDF5Constants.H5P_DEFAULT);

			// Close everything
			int status = H5.H5Fclose(fileId);

		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return hasData;
	}

	/**
	 * This operation reads the concentration at a specific position from the
	 * given HDF5 file
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which concentrations will be
	 *            copied
	 * @param lastTimeStep
	 *            The value of the last time step to know which concentration
	 *            group to look at
	 * @param i
	 *            The index for the position on the x axis
	 * @param j
	 *            The index for the position on the y axis
	 * @param k
	 *            The index for the position on the z axis
	 * @return A two dimensional array containing the concentrations
	 */
	public double[][] readConcentration(String fromName, int lastTimeStep, int i, int j, int k) {
		// The array for the concentrations
		double[][] concentration = new double[0][0];

		try {
			// Open the file
			int fileId = H5.H5Fopen(fromName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);

			// Create the name of the dataset
			String datasetName = "concentrationsGroup/concentration_" + lastTimeStep + "/position_" + i + "_" + j + "_"
					+ k;

			// Open it
			int datasetId = H5.H5Dopen(fileId, datasetName, HDF5Constants.H5P_DEFAULT);

			// Get the dataspace object
			int dataspaceId = H5.H5Dget_space(datasetId);

			// Get the dimensions of the dataset
			long[] dims = new long[2];
			int status = H5.H5Sget_simple_extent_dims(dataspaceId, dims, null);

			// Create the array that will receive the concentrations
			concentration = new double[(int) dims[0]][(int) dims[1]];

			// Read the data set
			status = H5.H5Dread(datasetId, HDF5Constants.H5T_IEEE_F64LE, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL,
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
	 * This operation writes the concentrations at this position in the new HDF5
	 * file
	 * 
	 * @param toName
	 *            The name of the created HDF5 file
	 * @param i
	 *            The index for the position on the x axis
	 * @param j
	 *            The index for the position on the y axis
	 * @param k
	 *            The index for the position on the z axis
	 * @param concentration
	 *            The two dimensional array containing the concentrations
	 * @param map
	 *            The map containing the index and composition of the clusters
	 *            in the previous network
	 * @param clusters
	 *            The list of clusters composing the new network
	 */
	public void writeConcentration(String toName, int i, int j, int k, double[][] concentration, int[][] map,
			ArrayList<Cluster> clusters) {
		// The list storing the new index and concentrations
		ArrayList<double[]> newConc = new ArrayList<double[]>();

		try {
			// Looping over the array on concentration to check that the
			// clusters
			// from the previous network still exist in the new network
			for (int l = 0; l < concentration.length; l++) {
				// boolean to know if the cluster is in the new network
				boolean found = false;
				// index in the new network
				int index = 0;
				// Loop on the new network
				for (Cluster cluster : clusters) {
					// Check the composition
					if ((cluster.nHe == map[(int) concentration[l][0]][0] || cluster.nXe == map[(int) concentration[l][0]][0])
							&& (cluster.nV == map[(int) concentration[l][0]][1])
							&& (cluster.nI == map[(int) concentration[l][0]][2])) {
						// Add the cluster to the new list, with the new index
						// but previous concentration
						double[] concCluster = { (double) index, concentration[l][1] };
						newConc.add(concCluster);
						found = true;
					}
					// increment the index
					index++;
				}

				// If the cluster was not found
				if (!found) {
					// Inform the user
					System.out.println("Cluster with the following composition (He/Xe, V, I): "
							+ map[(int) concentration[l][0]][0] + ", " + map[(int) concentration[l][0]][1] + ", "
							+ map[(int) concentration[l][0]][2]
							+ " is not present in the new network. Its concentration was " + concentration[l][1]
							+ " at the (" + i + ", " + j + ", " + k + ")-th grid point.");
				}
			}

			// Recreate the concentrattion array from the new concentration list
			concentration = new double[newConc.size()][2];
			// Loop on it to fill it
			for (int l = 0; l < concentration.length; l++) {
				double[] temp = newConc.get(l);
				concentration[l][0] = temp[0];
				concentration[l][1] = temp[1];
			}

			// Open the created HDF5 file
			int fileId = H5.H5Fopen(toName, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);

			// Create the same dataset in the new file
			String datasetName = "concentrationsGroup/concentration_0/position_" + i + "_" + j + "_" + k;

			// Create the dataspace for the dataset with dimension dims
			long[] dims = new long[2];
			dims[0] = concentration.length;
			dims[1] = 2;
			int concDataspaceId = H5.H5Screate_simple(2, dims, null);

			// Create the dataset for the concentrations
			int datasetId = H5.H5Dcreate(fileId, datasetName, HDF5Constants.H5T_IEEE_F64LE, concDataspaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Write the concentration array in the dataset
			int status = H5.H5Dwrite(datasetId, HDF5Constants.H5T_IEEE_F64LE, HDF5Constants.H5S_ALL,
					HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, concentration);

			// Close everything
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
	 * This operation copies the concentrations from a given HDF5 file to the
	 * one created by the preprocessor
	 * 
	 * @param fromName
	 *            The name of the HDF5 file from which the concentrations will
	 *            be copied
	 * @param toName
	 *            The name of the created HDF5 file
	 * @param clusters
	 *            The list of clusters that compose the network that will be
	 *            written in the created HDF5 file
	 */
	public void copyConcentration(String fromName, String toName, int[] gridSize, ArrayList<Cluster> clusters) {
		try {
			// Read the concentration group from the given file to get the value
			// of the last timestep
			int lastTimeStep = readLastTimeStep(fromName);

			if (lastTimeStep >= 0) {
				// Create the concentration group in the new HDF5 file
				createInitialConcentrationGroup(toName);

				// Read the times from the given HDF5 file
				double[] times = readTimes(fromName, lastTimeStep);

				// Create the concentration subgroup in the new file and store
				// the times
				createConcentrationSubGroup(toName, times);

				// Copy the surface position information
				copySurface(fromName, lastTimeStep, toName);

				// Create a map of the index and composition of the network
				// present
				// in the given HDF5 file
				int[][] networkMap = readNetwork(fromName);

				// Loop on all the position to read and copy the values of the
				// concentrations
				for (int k = -1; k <= gridSize[2]; k++) {
					for (int j = -1; j <= gridSize[1]; j++) {
						for (int i = 0; i <= gridSize[0]; i++) {
							if (!hasConcentrationDataset(fromName, lastTimeStep, i, j, k))
								continue;

							// Read the concentrations at this position from the
							// given HDF5 file
							double[][] conc = readConcentration(fromName, lastTimeStep, i, j, k);

							// Write the concentrations for this position in the
							// new HDF5 file
							writeConcentration(toName, i, j, k, conc, networkMap, clusters);
						}
					}
				}
			}

			else {
				// Complain
				System.err.println("The given HDF5 file doesn't have stored concentration. Aborting.");
				return;
			}

		} catch (RuntimeException e) {
			// Complain
			e.printStackTrace();
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
		try {
			// Open the HDF5 file
			int fileId = H5.H5Fopen(name, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT);

			// Create the header group
			int networkGroupId = H5.H5Gcreate(fileId, "networkGroup", HDF5Constants.H5P_DEFAULT,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Create the array that will store the network
			int networkSize = clusters.size();
			double[][] networkArray = new double[networkSize][6];

			int id = 0;
			// Loop on the clusters
			for (Cluster cluster : clusters) {
				// Store the composition
				networkArray[id][0] = cluster.nHe;
				if (maxXe > 0)
					networkArray[id][0] = cluster.nXe;
				networkArray[id][1] = cluster.nV;
				networkArray[id][2] = cluster.nI;

				// Store the formation energy
				networkArray[id][3] = cluster.E_f;

				// Store the migration energy
				networkArray[id][4] = cluster.E_m;

				// Store the diffusion factor
				networkArray[id][5] = cluster.D_0;

				// increment the id number
				id++;
			}

			// Create the dataspace for the network with dimension dims
			long[] dims = new long[2];
			dims[0] = networkSize;
			dims[1] = 6;
			int networkDataSpaceId = H5.H5Screate_simple(2, dims, null);

			// Create the dataset for the network
			int datasetId = H5.H5Dcreate(networkGroupId, "network", HDF5Constants.H5T_IEEE_F64LE, networkDataSpaceId,
					HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Write networkArray in the dataset
			int status = H5.H5Dwrite(datasetId, HDF5Constants.H5T_IEEE_F64LE, HDF5Constants.H5S_ALL,
					HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, networkArray);

			// Create the attribute for the network size
			int networkSizeDataSpaceId = H5.H5Screate(HDF5Constants.H5S_SCALAR);
			int networkSizeAttributeId = H5.H5Acreate(datasetId, "networkSize", HDF5Constants.H5T_STD_I32LE,
					networkSizeDataSpaceId, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);

			// Write it
			int[] tempNetworkSize = { networkSize };
			status = H5.H5Awrite(networkSizeAttributeId, HDF5Constants.H5T_STD_I32LE, tempNetworkSize);

			// Close everything
			status = H5.H5Aclose(networkSizeAttributeId);
			status = H5.H5Dclose(datasetId);
			status = H5.H5Gclose(networkGroupId);
			status = H5.H5Fclose(fileId);
		} catch (RuntimeException e) {
			// Complain
			e.printStackTrace();
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}

		return;
	}
}