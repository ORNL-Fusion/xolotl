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

		// Get command line arguments
		try {

			myArgs = CliFactory.parseArguments(Arguments.class, args);

			if (myArgs != null) {
				try {
					// Create the Preprocessor
					Preprocessor preprocessor = new Preprocessor(myArgs);

					// Generate the network of clusters
					ArrayList<Cluster> clusters = preprocessor
							.generateNetwork(args);

					// Create the HDF5 file
					preprocessor.createHDF5("networkInit.h5");

					// Write the header in it
					int[] dim = { 8 };
					int[] refinement = { 0 };
					preprocessor.writeHeader("networkInit.h5", dim, refinement);

					// Write the network in it
					preprocessor.writeNetwork("networkInit.h5", clusters);

					// Write the file containing the parameters that are needed
					// to run Xolotl
					preprocessor.writeParameterFile("params.txt",
							preprocessor.xolotlParams);

				} catch (Exception e) {
					e.printStackTrace();
				}
			}

		} catch (ArgumentValidationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		return;
	}

}
