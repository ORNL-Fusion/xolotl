package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import java.util.Enumeration;
import java.util.Properties;
import java.util.ArrayList;
import java.io.*;

import gov.ornl.xolotl.preprocessor.Preprocessor;
import gov.ornl.xolotl.preprocessor.Arguments;
import gov.ornl.xolotl.preprocessor.Cluster;

import org.junit.Test;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;

/**
 * This class is responsible for testing the Preprocessor class
 */
public class PreprocessorTest {

	/**
	 * This operation checks that the default parameters will be used along with
	 * writeParameterFile and loadParameterFile.
	 */
	@Test
	public void testParameterFile() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] {});

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Write the parameter file
				preprocessor.writeParameterFile("paramsTest",
						preprocessor.xolotlParams);

				// Load the properties from the parameter file to check they
				// were written correctly
				Properties inProps = preprocessor
						.loadParameterFile("paramsTest");
				// Enumeration to hold the parameter names
				Enumeration<?> paramNames = inProps.propertyNames();
				while (paramNames.hasMoreElements()) {
					String key = (String) paramNames.nextElement();
					String value = inProps.getProperty(key);
					// Check that the default parameter values were used
					assertEquals(preprocessor.xolotlParams.getProperty(key),
							value);
				}

				// Delete the parameter file
				new File("paramsTest").delete();

			}
		} catch (ArgumentValidationException e1) {
			e1.printStackTrace();
		}

		return;
	}

	/**
	 * This operation checks that the options specified via the command line
	 * will override the default values.
	 */
	@Test
	public void testCLOptionOverride() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--perfHandler", "dummy",
							"--petscArgs=" + "-da_grid_x 8 -ts_final_time 2" });

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Write the parameter file
				preprocessor.writeParameterFile("clOptionsTest",
						preprocessor.xolotlParams);

				// Load the properties from the parameter file to check they
				// were written correctly
				Properties inProps = preprocessor
						.loadParameterFile("clOptionsTest");
				// Enumeration to hold the parameter names
				Enumeration<?> paramNames = inProps.propertyNames();
				while (paramNames.hasMoreElements()) {
					String key = (String) paramNames.nextElement();
					String value = inProps.getProperty(key);
					// Check that the default parameter values were used
					assertEquals(preprocessor.xolotlParams.getProperty(key),
							value);
				}

				// Delete the parameter file
				new File("clOptionsTest").delete();

			}
		} catch (ArgumentValidationException e1) {
			e1.printStackTrace();
		}

		return;
	}

	/**
	 * This operation checks if the optional options are specified via the
	 * command line, that they will be included in the parameter file.
	 */
	@Test
	public void testOptionalOptions() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--material", "Fe", "--startTemp", "900" });

			// Check if there is a material argument
			assertEquals(true, parsedArgs.isMaterial());
			
			// Check that the material is Fe
			assertEquals("Fe", parsedArgs.getMaterial());
			
			// Check that the startTemp is 900
			assertEquals("900", parsedArgs.getStartTemp());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Write the parameter file
				preprocessor.writeParameterFile("optionalOpsTest",
						preprocessor.xolotlParams);

				// Load the properties from the parameter file to check they
				// were written correctly
				Properties inProps = preprocessor
						.loadParameterFile("optionalOpsTest");
				// Enumeration to hold the parameter names
				Enumeration<?> paramNames = inProps.propertyNames();
				while (paramNames.hasMoreElements()) {
					String key = (String) paramNames.nextElement();
					String value = inProps.getProperty(key);
					// Check that the default parameter values were used
					assertEquals(preprocessor.xolotlParams.getProperty(key),
							value);
				}

				// Delete the parameter file
				new File("optionalOpsTest").delete();

			}
		} catch (ArgumentValidationException e1) {
			e1.printStackTrace();
		}

		return;
	}

	/**
	 * This operation checks generateGrid.
	 */
	@Test
	public void testGenerateGrid() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] {});

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Create the known grid array
				double[] knownGrid = new double[15];

				// Fill it with known values
				for (int i = 0; i < 15; i++) {
					knownGrid[i] = (double) i * 0.57142857142;
				}

				try {
					// Generate a grid
					double[] newGrid = preprocessor.generateGrid(8, 1);

					// Check the length of it
					assertEquals(newGrid.length, knownGrid.length);

					// Check all the values
					for (int i = 0; i < newGrid.length; i++) {
						assertEquals(newGrid[i], knownGrid[i], 1.0e-5);
					}

				} catch (Exception e) {
					// Complain and fail
					e.printStackTrace();
					fail();
				}

			}
		} catch (ArgumentValidationException e1) {
			e1.printStackTrace();
		}

		return;
	}

	/**
	 * This operation checks the writing of the HDF5 file.
	 */
	@Test
	public void testHDF5Writing() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] {});

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Create an empty cluster array
				ArrayList<Cluster> clusters = new ArrayList<Cluster>();

				// Create a cluster
				Cluster cluster = new Cluster();
				cluster.nHe = 1;
				cluster.nV = 23;
				cluster.nI = 52;
				cluster.E_He = 9.325;
				cluster.E_V = 34.2346;
				cluster.E_I = 3326424.2323543;
				cluster.E_m = 0.04;
				cluster.D_0 = 1.1;

				// Add it to clusters
				clusters.add(cluster);

				try {
					// Create the HDF5 file
					preprocessor.createHDF5("test.h5");

					// Write the header in it
					int[] dim = { 4 };
					int[] refinement = { 0 };
					preprocessor.writeHeader("test.h5", dim, refinement);

					// Write the network in it
					preprocessor.writeNetwork("test.h5", clusters);

					// Check that the file was created
					File f = new File("test.h5");
					boolean fileExists = (f.exists() && !f.isDirectory());
					assertEquals(fileExists, true);

				} catch (Exception e) {
					// Complain and fail
					e.printStackTrace();
					fail();
				}

			}
		} catch (ArgumentValidationException e1) {
			e1.printStackTrace();
		}

		return;
	}
}
