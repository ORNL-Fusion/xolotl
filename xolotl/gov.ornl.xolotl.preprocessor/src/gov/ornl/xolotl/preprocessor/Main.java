/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.util.ArrayList;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;

/**
 * This class launches the Xolotl preprocessor.
 * 
 * @author Jay Jay Billings
 * 
 */
public class Main {

	/**
	 * This operation launches the preprocessor and creates the initial
	 * conditions for Xolotl.
	 * 
	 * @param args
	 *            Command line arguments.
	 */
	public static void main(String[] args) {
		// Local Declarations
		Arguments myArgs = null;

		// Print some information to let the user know what is happening.
		System.out.println("Starting Xolotl Preprocessor...");

		// Get command line arguments
		try {
			// Parse the command line arguments.
			myArgs = CliFactory.parseArguments(Arguments.class, args);
			System.out.println("Command line arguments loaded.");

			// Generate the input if the arguments are valid.
			if (myArgs != null) {
				// Create the Preprocessor
				Preprocessor preprocessor = new Preprocessor(myArgs);

				// Generate the network of clusters
				ArrayList<Cluster> clusters = preprocessor
						.generateNetwork(args);
				System.out.println("Network generated.");
				
				// Get the name of the networkFile from xolotlParams
				String networkFileName = preprocessor.xolotlParams.getProperty("networkFile");

				// Create the HDF5 file
				preprocessor.createHDF5(networkFileName);

				// Write the network in it
				preprocessor.writeNetwork(networkFileName, clusters);

				if (myArgs.isCheckpoint()) {
					String HDF5FileName = myArgs.getCheckpoint();
					// Read the header and the concentration from this file 
					// and copy them to the network file
					int gridPoints = preprocessor.copyHeader(HDF5FileName, networkFileName);
					preprocessor.copyConcentration(HDF5FileName, networkFileName, gridPoints);
				}
				else {
					// Write the header in it
					int[] dim = { 8 };
					int[] refinement = { 0 };
					preprocessor.writeHeader(networkFileName, dim, refinement);
				}

				System.out.println("HDF5 file generated.");

				// Write the file containing the parameters that are needed
				// to run Xolotl
				preprocessor.writeParameterFile("params.txt",
						preprocessor.xolotlParams);
				System.out.println("Parameters written.");
			}
		} catch (ArgumentValidationException e1) {
			System.err.println("Unable to print input arguments. Aborting.");
			e1.printStackTrace();
			return;
		} catch (Exception e) {
			System.err.println("Exception caught while generating input. "
					+ "Aborting.");
			e.printStackTrace();
			return;
		}

		// Say goodbye
		System.out.println("Finished.");

		return;
	}
}
