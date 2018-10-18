package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import java.util.Enumeration;
import java.util.Properties;
import java.util.ArrayList;
import java.io.*;

import gov.ornl.xolotl.preprocessor.Preprocessor;
import gov.ornl.xolotl.preprocessor.Arguments;

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
			parsedArgs = CliFactory.parseArguments(Arguments.class, new String[] {});

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Write the parameter file
				preprocessor.writeParameterFile("paramsTest", preprocessor.xolotlParams);

				// Load the properties from the parameter file to check they
				// were written correctly
				Properties inProps = preprocessor.loadParameterFile("paramsTest");

				// Enumeration to hold the parameter names
				Enumeration<?> paramNames = inProps.propertyNames();
				while (paramNames.hasMoreElements()) {
					String key = (String) paramNames.nextElement();
					String value = inProps.getProperty(key);
					// Check that the default parameter values were used
					assertEquals(preprocessor.xolotlParams.getProperty(key), value);
				}

				// Delete the parameter file
				new File("paramsTest").delete();
			}
		} catch (ArgumentValidationException e) {
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
					new String[] { "--perfHandler", "dummy", "--petscArgs=" + "-da_grid_x 8 -ts_final_time 2" });

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Write the parameter file
				preprocessor.writeParameterFile("clOptionsTest", preprocessor.xolotlParams);

				// Load the properties from the parameter file to check they
				// were written correctly
				Properties inProps = preprocessor.loadParameterFile("clOptionsTest");

				// Enumeration to hold the parameter names
				Enumeration<?> paramNames = inProps.propertyNames();
				while (paramNames.hasMoreElements()) {
					String key = (String) paramNames.nextElement();
					String value = inProps.getProperty(key);
					// Check that the default parameter values were used
					assertEquals(preprocessor.xolotlParams.getProperty(key), value);
				}

				// Delete the parameter file
				new File("clOptionsTest").delete();
			}
		} catch (ArgumentValidationException e) {
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
					new String[] { "--startTemp", "900", "--flux", "1.5" });

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);

				// Write the parameter file
				preprocessor.writeParameterFile("optionalOpsTest", preprocessor.xolotlParams);

				// Load the properties from the parameter file to check they
				// were written correctly
				Properties inProps = preprocessor.loadParameterFile("optionalOpsTest");

				// Enumeration to hold the parameter names
				Enumeration<?> paramNames = inProps.propertyNames();
				while (paramNames.hasMoreElements()) {
					String key = (String) paramNames.nextElement();
					String value = inProps.getProperty(key);
					// Check that the default parameter values were used
					assertEquals(preprocessor.xolotlParams.getProperty(key), value);
				}

				// Delete the parameter file
				new File("optionalOpsTest").delete();
			}
		} catch (ArgumentValidationException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		}

		return;
	}

	/**
	 * This operation checks that it is not possible to give wrong sizes for He,
	 * I, and V.
	 */
	@Test(expected = IllegalArgumentException.class)
	public void testBadMaxClusterSizeOptions() {

		// Local Declarations
		Arguments parsedArgs = null;

		try {
			// Try a number of helium that is too big
			parsedArgs = CliFactory.parseArguments(Arguments.class, new String[] { "--maxHeSize", "10" });

			// Check that the max helium cluster size is 10
			assertEquals(10, parsedArgs.getMaxHeSize());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);
				fail("Should have thrown an IllegalArgumentException because "
						+ "the maximum He size is out of range.");
			}
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}

		try {
			// Try a number of vacancy that is negative
			parsedArgs = CliFactory.parseArguments(Arguments.class, new String[] { "--maxVSize", "-1" });

			// Check that the max interstitial cluster size is -1
			assertEquals(-2, parsedArgs.getMaxVSize());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);
				fail("Should have thrown an IllegalArgumentException because " + "the maximum V size is out of range.");
			}
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}

		try {
			// Try a number of interstitial that is negative
			parsedArgs = CliFactory.parseArguments(Arguments.class, new String[] { "--maxISize", "-2" });

			// Check that the max interstitial cluster size is -2
			assertEquals(-2, parsedArgs.getMaxISize());

			if (parsedArgs != null) {
				Preprocessor preprocessor = new Preprocessor(parsedArgs);
				fail("Should have thrown an IllegalArgumentException because " + "the maximum I size is out of range.");
			}
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}

		return;
	}
}