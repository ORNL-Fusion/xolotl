/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.util.ArrayList;

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

		// Create the Preprocessor
		Preprocessor preprocessor = new Preprocessor();
		
		// Generate the clusters
		ArrayList<Cluster> clusters = preprocessor.generate(args);
		
		// Dump the clusters to stdout
		for (Cluster cluster : clusters) {
			System.out.println(cluster.toString());
		}
		
		return;
	}

}
