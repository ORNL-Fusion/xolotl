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

				// Create the HDF5 file
				preprocessor.createHDF5("networkInit.h5");

				// Write the header in it
				int[] dim = { 8 };
				int[] refinement = { 0 };
				preprocessor.writeHeader("networkInit.h5", dim, refinement);

				// Write the network in it
				preprocessor.writeNetwork("networkInit.h5", clusters);
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
