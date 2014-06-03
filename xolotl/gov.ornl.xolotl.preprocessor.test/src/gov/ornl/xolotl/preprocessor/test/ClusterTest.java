/**
 * 
 */
package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;
import gov.ornl.xolotl.preprocessor.Cluster;

import org.junit.Test;

/**
 * This class is responsible for testing the Cluster class, specifically its
 * toString operation.
 * 
 * @author Jay Jay Billings
 * 
 */
public class ClusterTest {

	/**
	 * Test method for {@link gov.ornl.xolotl.preprocessor.Cluster#toString()}.
	 */
	@Test
	public final void testToString() {

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

		// Create the test string and make sure they match.
		String clusterAsString = "1 23 52 9.325 34.2346 3326424.2323543 0.04 1.1";
		assertEquals(clusterAsString, cluster.toString());

		// Modify the binding energies to be infinite and check the strings
		// again. The cluster class should return "Infinity" instead of a value
		// for the energies.
		cluster.E_He = Double.POSITIVE_INFINITY;
		cluster.E_V = Double.POSITIVE_INFINITY;
		cluster.E_I = Double.POSITIVE_INFINITY;
		cluster.E_m = Double.POSITIVE_INFINITY;
		clusterAsString = "1 23 52 Infinity Infinity Infinity Infinity 1.1";
		assertEquals(clusterAsString, cluster.toString());
		// Consequently, such a cluster would be the ultimate sink in the
		// network!

		// Modify the migration energy and diffusion factor and check that this
		// does not interfere with the other properties. 
		cluster.E_He = Double.POSITIVE_INFINITY;
		cluster.E_V = Double.POSITIVE_INFINITY;
		cluster.E_I = Double.POSITIVE_INFINITY;
		cluster.E_m = 0.13;
		cluster.D_0 = 2.9e-8;
		clusterAsString = "1 23 52 Infinity Infinity Infinity 0.13 2.9E-8";
		assertEquals(clusterAsString, cluster.toString());

		return;
	}

}
