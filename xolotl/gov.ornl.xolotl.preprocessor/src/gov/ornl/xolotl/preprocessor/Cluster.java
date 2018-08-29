/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

/**
 * A simple, all public constructor that represents a cluster. It overwrites
 * toString() to write itself according to the specification needed by Xolotl.
 * 
 * @author jaybilly
 * 
 */
public class Cluster {

	// The number of helium atoms in the cluster. Default value is 0.
	public int nHe = 0;

	// The number of xenon atoms in the cluster. Default value is 0.
	public int nXe = 0;

	// The number of vacancies in the cluster. Default value is 0.
	public int nV = 0;

	// The number of interstitials in the cluster. Default value is 0.
	public int nI = 0;

	// The number of deuterium atoms in the cluster. Default value is 0.
	public int nD = 0;

	// The number of tritium atoms in the cluster. Default value is 0.
	public int nT = 0;

	// The formation energy of this cluster
	public double E_f = 0.0;

	// The migration energy with which the cluster diffuses. Default value is
	// infinity so that the cluster is immobile.
	public double E_m = Double.POSITIVE_INFINITY;

	// The diffusion factor that sets the scale of the diffusion. Default value
	// is infinity so that the cluster is immobile.
	// is 0.
	public double D_0 = 0.0;

}
