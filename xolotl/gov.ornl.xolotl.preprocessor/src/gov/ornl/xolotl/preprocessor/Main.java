/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.io.File;
import java.util.ArrayList;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;
import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.h5.H5File;

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

		// Get command line arguments
		try {
			Arugments myArgs = CliFactory.parseArguments(Arugments.class, args);
		} catch (ArgumentValidationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// Create a file from the uri
		File file = new File("test.hdf5");

		// Retrieve an instance of the HDF5 format
		FileFormat fileFormat = FileFormat
				.getFileFormat(FileFormat.FILE_TYPE_HDF5);

		// Create an H5 file. If it exists already, then delete it.
		try {
			H5File h5File = (H5File) fileFormat.createFile(file.getPath(),
					FileFormat.FILE_CREATE_DELETE);
		} catch (Exception e) {
			// Complain
			e.printStackTrace();
		}
		
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
