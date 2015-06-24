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
import org.junit.rules.ExpectedException;

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
				Properties inProps = preprocessor.loadParameterFile("paramsTest");
				
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
		}
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
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
				Properties inProps = preprocessor.loadParameterFile("clOptionsTest");
			
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
		}
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
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
					new String[] { "--startTemp", "900", "--heFlux", "1.5" });

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
		}
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		} 

		return;
	}
	
	/**
	 * This operation checks that it is not possible to give wrong 
	 * sizes for He, I, and V.
	 */
	@Test(expected=IllegalArgumentException.class)
	public void testBadMaxClusterSizeOptions() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			// Try a number of helium that is too big
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--maxHeSize", "10" });
			
			// Check that the max helium cluster size is 10
			assertEquals(10, parsedArgs.getMaxHeSize());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);
				fail("Should have thrown an IllegalArgumentException because "
						+ "the maximum He size is out of range.");
			}
		}
		catch (ArgumentValidationException e) {
			e.printStackTrace();
		} 

		try {
			// Try a number of vacancy that is negative
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--maxVSize", "-1" });
			
			// Check that the max interstitial cluster size is -1
			assertEquals(-2, parsedArgs.getMaxVSize());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);
				fail("Should have thrown an IllegalArgumentException because "
						+ "the maximum V size is out of range.");
			}
		}
		catch (ArgumentValidationException e) {
			e.printStackTrace();
		} 
		
		try {
			// Try a number of interstitial that is negative
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--maxISize", "-2" });
			
			// Check that the max interstitial cluster size is -2
			assertEquals(-2, parsedArgs.getMaxISize());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);
				fail("Should have thrown an IllegalArgumentException because "
						+ "the maximum I size is out of range.");
				}
		}
		catch (ArgumentValidationException e) {
			e.printStackTrace();
		} 

		return;
	}

	/**
	 * This operation checks the generation of the network.
	 */
	@Test
	public void testNetworkGeneration() {
		// Local Declarations
		Arguments parsedArgs = null;

		try {
			// Keep the default arguments
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] {});

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Generate the network
				ArrayList<Cluster> network = preprocessor.generateNetwork();

				// Check the size of the network
				assertEquals(network.size(), 2067);
			}

			// Change the number of V clusters
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--maxVSize", "60" });

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Generate the network
				ArrayList<Cluster> network = preprocessor.generateNetwork();

				// Check the size of the network
				assertEquals(network.size(), 7678);
			}

			// Use the phase-cut method
			parsedArgs = CliFactory.parseArguments(Arguments.class,
					new String[] { "--phaseCut"});

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Generate the network
				ArrayList<Cluster> network = preprocessor.generateNetwork();

				// Check the size of the network
				assertEquals(network.size(), 476);
			}
		}
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
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
				cluster.E_f = 12.3;
				cluster.E_m = 0.04;
				cluster.D_0 = 1.1;

				// Add it to clusters
				clusters.add(cluster);

				// Create the HDF5 file
				preprocessor.createHDF5("test.h5");

				// Write the header in it
				preprocessor.writeHeader("test.h5", parsedArgs);

				// Write the network in it
				preprocessor.writeNetwork("test.h5", clusters);

				// Check that the file was created
				File f = new File("test.h5");
				boolean fileExists = (f.exists() && !f.isDirectory());
				assertEquals(fileExists, true);
			}
		}
		catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		} 

		return;
	}
}