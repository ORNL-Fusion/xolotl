package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import gov.ornl.xolotl.preprocessor.FormationEnergyEngine;

import org.junit.Test;

/**
 * This class is responsible for testing the FormationEnergyEngine class to make
 * sure that it correctly reports formation energies for clusters.
 * 
 * @author Jay Jay Billings
 * 
 */
public class FormationEnergyEngineTest {

	/**
	 * This operation checks the Helium formation energies.
	 */
	@Test
	public void testHeEnergies() {

		// Local Declarations
		FormationEnergyEngine engine = new FormationEnergyEngine();

		// Check the code for He_2.
		assertEquals(11.44, engine.getHeFormationEnergy(2), 1.0e-2);

		// Check the code for He_8.
		assertEquals(38.80, engine.getHeFormationEnergy(8), 1.0e-2);

		// Check the code for the invalid case with size = 0.
		assertEquals(Double.POSITIVE_INFINITY, engine.getHeFormationEnergy(0),
				1.0e-16);

		// Check the code for the invalid case where size > 8.
		assertEquals(Double.POSITIVE_INFINITY, engine.getHeFormationEnergy(10),
				1.0e-16);

		return;
	}

	/**
	 * This operation checks the vacancy formation energies.
	 */
	@Test
	public void testVEnergies() {

		// Local Declarations
		FormationEnergyEngine engine = new FormationEnergyEngine();

		// Check the code for V_2.
		assertEquals(7.25, engine.getVFormationEnergy(2), 1.0e-2);

		// Check the code for the invalid case with size = 0.
		assertEquals(Double.POSITIVE_INFINITY, engine.getVFormationEnergy(0),
				1.0e-16);

		// Check V = 29 to make sure the larger size clusters are computed
		// correctly.
		assertEquals(51.765, engine.getVFormationEnergy(29), 1.0e-2);

		return;
	}

	/**
	 * This operation checks the interstitial formation energies.
	 */
	@Test
	public void testIEnergies() {

		// Local Declarations
		FormationEnergyEngine engine = new FormationEnergyEngine();

		// Check the code for I_2.
		assertEquals(18.5, engine.getIFormationEnergy(2), 1.0e-2);

		// Check the code for I_6.
		assertEquals(48.0, engine.getIFormationEnergy(6), 1.0e-2);

		// Check the code for I_10.
		assertEquals(72.0, engine.getIFormationEnergy(10), 1.0e-16);

		// Check the code for the invalid case with size = 0.
		assertEquals(Double.POSITIVE_INFINITY, engine.getIFormationEnergy(0),
				1.0e-16);

		return;
	}

	/**
	 * This operation checks the Helium-Vacancy formation energies.
	 * 
	 * The fit coefficients and sample values of the function are available in
	 * the Xolotl repository. They were used to compute the test values.
	 * 
	 */
	@Test
	public void testHeVEnergies() {
		// Local Declarations
		FormationEnergyEngine engine = new FormationEnergyEngine();

		// Check the code for He_30V_6
		assertEquals(75.049, engine.getHeVFormationEnergy(30, 6), 1.0e-2);

		// Check the code for He_1V_44
		assertEquals(66.311, engine.getHeVFormationEnergy(1, 44), 1.0e-2);

		// Check the code for He_1V_1
		assertEquals(5.142, engine.getHeVFormationEnergy(1, 1), 1.0e-2);

		return;
	}

	/**
	 * This operation makes sure that the FormationEnergyEngine can read
	 * coefficients from a file during construction.
	 */
	@Test
	public void testFitFile() {

		try {
			// Create the fit file
			File fitFile = new File("fit.csv");
			FileWriter fitFileWriter = new FileWriter(fitFile);
			BufferedWriter fitWriter = new BufferedWriter(fitFileWriter);

			// Write to the file
			for (int i = 0; i < 8; i++) {
				fitWriter.write("0.0,0.0,0.0,0.0,0.0,0.0\n");
			}
			fitWriter.close();

			// Create the FormationEnergyEngine and test it. The answer should be
			// equal to 0.0.
			FormationEnergyEngine engine = new FormationEnergyEngine();
			assertEquals(0.0, engine.getHeVFormationEnergy(31, 6), 1.0e-2);

			// Delete the fit file
			fitFile.delete();

			// Check the code for He_30V_6 -> He_29V_6 + He_1 using the original
			// coefficients. (Make sure the values from the file don't hang
			// around.)
			engine = new FormationEnergyEngine();
			assertEquals(75.049, engine.getHeVFormationEnergy(30, 6), 1.0e-2);

		} catch (IOException e) {
			// Complain and fail
			e.printStackTrace();
			fail();
		}

		return;
	}

}
