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

/**
 * This class generates a valid parameter file for Xolotl.
 * 
 */
public class Preprocessor {

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
	 * This operation generates a single string containing the PETSc arguments that
	 * will be used to set the PETSc parameter required to run Xolotl.
	 * 
	 * @param petscArgs The PETSc command line arguments as specified by the
	 *                  Arguments interface
	 * @return A string containing the PETSc arguments that will be passed to Xolotl
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
		petscOptions.put("-ts_adapt_wnormtype", "INFINITY");

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
	 * @param args The command line arguments as specified by the Arguments
	 *             interface.
	 */
	public Preprocessor(Arguments args) {
		// Set the maximum size of a helium cluster in the network.
		int maxHe = args.getMaxHeSize();
		// Check to make sure the user entered an appropriate value
		if ((maxHe > 8) || (maxHe < 0)) {
			throw new IllegalArgumentException("The maxium helium size must be less than 9 ( 0 <= maxHeSize < 9 )");
		}

		// Set the maximum size of a xenon cluster in the network.
		int maxXe = args.getMaxXeSize();
		if (maxXe < 0) {
			throw new IllegalArgumentException("The maxium xenon size must be positive ( 0 <= maxXeSize )");
		}

		// Set the maximum size of a deuterium cluster in the network.
		int maxD = args.getMaxDSize();
		if (maxD < 0) {
			throw new IllegalArgumentException("The maxium deuterium size must be positive ( 0 <= maxDSize )");
		}

		// Set the maximum size of a tritium cluster in the network.
		int maxT = args.getMaxTSize();
		if (maxT < 0) {
			throw new IllegalArgumentException("The maxium tritium size must be positive ( 0 <= maxTSize )");
		}

		// Set the maximum size of a vacancy cluster in the network.
		int maxV = args.getMaxVSize();
		if (maxV < 0) {
			throw new IllegalArgumentException("The maxium vacancy must be positive ( 0 <= maxVSize )");
		}

		// The maximum size of an interstitial cluster in the network.
		int maxI = args.getMaxISize();
		if (maxI < 0) {
			throw new IllegalArgumentException("The maxium interstitial must be positive ( 0 <= maxISize )");
		}

		// The number of dimension for the problem to solve
		dim = Integer.parseInt(args.getDimensions());

		// Set the parameter options that will be passed to Xolotl
		xolotlParams.setProperty("dimensions", args.getDimensions());
		xolotlParams.setProperty("startTemp", args.getStartTemp());
		xolotlParams.setProperty("flux", args.getFlux());
		xolotlParams.setProperty("material", args.getMaterial());
		xolotlParams.setProperty("perfHandler", args.getPerfHandler());
		xolotlParams.setProperty("vizHandler", args.getVizHandler());
		xolotlParams.setProperty("petscArgs", generatePetscArgs(args.getPetscArgs()));
		xolotlParams.setProperty("process", args.getProcess());

		// The following parameter options are optional and will only
		// be set if they are specified via the command line
		if (args.isTempFile())
			xolotlParams.setProperty("tempFile", args.getTempFile());
		if (args.isHeat())
			xolotlParams.setProperty("heat", args.getHeat());
		if (args.isFluxFile())
			xolotlParams.setProperty("fluxFile", args.getFluxFile());
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
		if (args.isBurstingDepth())
			xolotlParams.setProperty("burstingDepth", args.getBurstingDepth());

		// Build the network argument
		String netString;
		if (maxXe > 0)
			netString = Integer.toString(args.getMaxXeSize());
		else
			netString = Integer.toString(args.getMaxHeSize()) + " " + Integer.toString(args.getMaxDSize()) + " "
					+ Integer.toString(args.getMaxTSize()) + " " + Integer.toString(args.getMaxVSize()) + " "
					+ Integer.toString(args.getMaxISize()) + " " + Boolean.toString(args.isPhaseCut());
		xolotlParams.setProperty("netParam", netString);

		// Build the grid argument
		String gridString = Integer.toString(args.getNxGrid()) + " " + Double.toString(args.getXStepSize()) + " "
				+ Integer.toString(args.getNyGrid()) + " " + Double.toString(args.getYStepSize()) + " "
				+ Integer.toString(args.getNzGrid()) + " " + Double.toString(args.getZStepSize());
		xolotlParams.setProperty("grid", gridString);

		// Set the default boundary conditions
		String boundaryString;
		if (maxXe > 0)
			boundaryString = "0 0";
		else
			boundaryString = "1 0";
		xolotlParams.setProperty("boundary", boundaryString);

	}

	/**
	 * This operation writes the parameters file that is needed to run Xolotl.
	 * 
	 * @param parameterFile The parameter file name
	 * @param parameters    The parameters that will be written to the file
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
	 * @param parameterFile The parameter file name
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
}