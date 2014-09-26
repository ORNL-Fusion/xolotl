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
		cluster.E_f = 9.325;
		cluster.E_m = 0.04;
		cluster.D_0 = 1.1;

		// Create the test string and make sure they match.
		String clusterAsString = "1 23 52 9.325 0.04 1.1";
		assertEquals(clusterAsString, cluster.toString());

		// Modify the formation energy to be infinite and check the strings
		// again. The cluster class should return "Infinity" instead of a value
		// for the energies.
		cluster.E_f = Double.POSITIVE_INFINITY;
		cluster.E_m = Double.POSITIVE_INFINITY;
		clusterAsString = "1 23 52 Infinity Infinity 1.1";
		assertEquals(clusterAsString, cluster.toString());
		// Consequently, such a cluster would be the ultimate sink in the
		// network!

		// Modify the migration energy and diffusion factor and check that this
		// does not interfere with the other properties. 
		cluster.E_f = Double.POSITIVE_INFINITY;
		cluster.E_m = 0.13;
		cluster.D_0 = 2.9e-8;
		clusterAsString = "1 23 52 Infinity 0.13 2.9E-8";
		assertEquals(clusterAsString, cluster.toString());

		return;
	}

}
