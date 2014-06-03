package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;
import gov.ornl.xolotl.preprocessor.BindingEnergyEngine;

import org.junit.Test;

/**
 * This class is responsible for testing the BindingEnergyEngine class to make
 * sure that it correctly reports binding energies for clusters.
 * 
 * @author Jay Jay Billings
 * 
 */
public class BindingEnergyEngineTest {

	/**
	 * This operation checks the Helium binding energies.
	 */
	@Test
	public void testHeEnergies() {

		// Local Declarations
		BindingEnergyEngine engine = new BindingEnergyEngine();

		// Check the code for He_2 --> He_1 + He_1.
		assertEquals(0.86, engine.getHeBindingEnergy(2), 1.0e-2);

		// Check the code for He_8 --> He_7 + He_1.
		assertEquals(2.28, engine.getHeBindingEnergy(8), 1.0e-2);

		// Check the code for the invalid case with size = 1.
		assertEquals(Double.POSITIVE_INFINITY, engine.getHeBindingEnergy(1),
				1.0e-16);

		// Check the code for the invalid case where size > 8.
		assertEquals(Double.POSITIVE_INFINITY, engine.getHeBindingEnergy(9),
				1.0e-16);

		return;
	}

	/**
	 * This operation checks the vacancy binding energies.
	 */
	@Test
	public void testVEnergies() {

		// Local Declarations
		BindingEnergyEngine engine = new BindingEnergyEngine();

		// Check the code for V_2 --> V_1 + V_1. That is the only valid size for
		// this cluster type at the moment.
		assertEquals(0.5519, engine.getVBindingEnergy(2), 1.0e-2);

		// Check the code for the invalid case with size = 1.
		assertEquals(Double.POSITIVE_INFINITY, engine.getVBindingEnergy(1),
				1.0e-16);

		// Check V = 29 to make sure the larger size clusters are computed
		// correctly.
		assertEquals(2.63,engine.getVBindingEnergy(29),1.0e-2);

		return;
	}

	/**
	 * This operation checks the interstitial binding energies.
	 */
	@Test
	public void testIEnergies() {

		// Local Declarations
		BindingEnergyEngine engine = new BindingEnergyEngine();

		// Check the code for I_2 --> I_1 + I_1.
		assertEquals(1.5, engine.getIBindingEnergy(2), 1.0e-2);

		// Check the code for I_6 --> I_5 + I_1.
		assertEquals(15.5, engine.getIBindingEnergy(6), 1.0e-2);

		// Check the code for I_10 --> I_9 + I_1.
		assertEquals(16.0, engine.getIBindingEnergy(10), 1.0e-16);

		// Check the code for the invalid case with size = 1.
		assertEquals(Double.POSITIVE_INFINITY, engine.getIBindingEnergy(1),
				1.0e-16);

		return;
	}

	/**
	 * This operation checks the Helium-Vacancy binding energies to single
	 * Helium and single vacancy clusters.
	 * 
	 * The fit coefficients and sample values of the function are available in
	 * the Xolotl repository. They were used to compute the test values.
	 * 
	 */
	@Test
	public void testHeVEnergies() {

		// Local Declarations
		BindingEnergyEngine engine = new BindingEnergyEngine();

		// Check the code for He_31V_6 -> He_30V_6 + He_1
		assertEquals(3.85, engine.getHeVtoHeBindingEnergy(31, 6), 1.0e-2);

		// Check the code for He_1V_44 -> He_1V_43 + V_1
		assertEquals(2.75, engine.getHeVtoVBindingEnergy(1, 44), 1.0e-2);

		// Check the code for He_8 --> He_7 + He_1 to make sure work is
		// delegated to the
		assertEquals(2.28, engine.getHeBindingEnergy(8), 1.0e-2);

		// Check that requesting He_1V_1 returns a value other than for this
		// special case.
		assertEquals(4.83, engine.getHeVtoHeBindingEnergy(1, 1), 1.0e-2);
		assertEquals(4.83, engine.getHeVtoVBindingEnergy(1, 1), 1.0e-2);

		// Check that larger V = 1 clusters like He_6V_1 are correct.
		assertEquals(2.56, engine.getHeVtoHeBindingEnergy(6, 1), 1.0e-2);
		assertEquals(11.78, engine.getHeVtoVBindingEnergy(6, 1), 1.0e-2);
		
		// Make sure that if there is no He or V the engine returns infinity for
		// the respective calls.
		assertEquals(Double.POSITIVE_INFINITY,
				engine.getHeVtoHeBindingEnergy(0, 1), 1.0e-2);
		assertEquals(Double.POSITIVE_INFINITY,
				engine.getHeVtoVBindingEnergy(1, 0), 1.0e-2);
		
		return;
	}

}
