package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import org.junit.Test;

import gov.ornl.xolotl.preprocessor.Arguments;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;

/**
 * This class is responsible for testing the Arguments command line interface.
 */
public class ArgumentsTest {

	/**
	 * This operation checks that the default argument values are used if no
	 * command line options are specified.
	 */
	@Test
	public void testDefaultArguments() {

		// Local Declarations
		final Arguments args;

		try {
			// Parse the empty string of arguments
			args = CliFactory.parseArguments(Arguments.class, new String[] {});

			// Check if there is a material argument
			assertEquals(false, args.isMaterial());

			// Check if there is a startTemp argument
			assertEquals("1000", args.getStartTemp());

			// Check if there is a tempFile argument
			assertEquals(false, args.isTempFile());

			// Check if there is an heFlux argument
			assertEquals(false, args.isHeFlux());
			
			// Check if there is an heFluence argument
			assertEquals(false, args.isHeFluence());

			// Check that the default perfHandler is std
			assertEquals("std", args.getPerfHandler());

			// Check if there is a vizHandler argument
			assertEquals("dummy", args.getVizHandler());

			// Check if there is a checkpoint argument
			assertEquals(false, args.isCheckpoint());

			// Check that the default networkFile is networkInit.h5
			assertEquals("networkInit.h5", args.getNetworkFile());

			// Check the default petscArgs
			assertEquals(
					"-da_grid_x 10 -ts_final_time 1000 "
							+ "-ts_max_steps 3 -ts_adapt_dt_max 10 -ts_max_snes_failures 200 "
							+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant "
							+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
					args.getPetscArgs());
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}
	}

	/**
	 * This operation tests that default parameter values are only overridden if
	 * they are specified via the command line and that the optional arguments 
	 * are only set if they are also specified
	 */
	@Test
	public void testSpecifiedArguments() {

		// Local Declarations
		final Arguments args;

		try {
			// Parse the specified string of arguments
			args = CliFactory.parseArguments(Arguments.class, new String[] {
					"--startTemp", "900", "--material", "Fe", "--perfHandler",
					"dummy" });
			
			// Check if there is a material argument
			assertEquals(true, args.isMaterial());
			
			// Check that the material is Fe
			assertEquals("Fe", args.getMaterial());
			
			// Check that the startTemp is 900
			assertEquals("900", args.getStartTemp());

			// Check if there is a tempFile argument
			assertEquals(false, args.isTempFile());

			// Check if there is an heFlux argument
			assertEquals(false, args.isHeFlux());
			
			// Check if there is an heFluence argument
			assertEquals(false, args.isHeFluence());

			// Check that the perfHandler is dummy
			assertEquals("dummy", args.getPerfHandler());

			// Check if there is a vizHandler argument
			assertEquals("dummy", args.getVizHandler());

			// Check if there is a checkpoint argument
			assertEquals(false, args.isCheckpoint());

			// Check that the default networkFile is networkInit.h5
			assertEquals("networkInit.h5", args.getNetworkFile());

			// Check the default petscArgs
			assertEquals(
					"-da_grid_x 10 -ts_final_time 1000 "
							+ "-ts_max_steps 3 -ts_adapt_dt_max 10 -ts_max_snes_failures 200 "
							+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant "
							+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
					args.getPetscArgs());
		} catch (ArgumentValidationException e) {
			e.printStackTrace();
		}
	}

}
