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

	// The number of Helium atoms in the cluster. Default value is 0.
	public int nHe = 0;

	// The number of vacancies in the cluster. Default value is 0.
	public int nV = 0;

	// The number of interstitials in the cluster. Default value is 0.
	public int nI = 0;

	// The binding energy of the cluster with He_1. Default value is 0.
	public double E_He = 0.0;

	// The binding energy of the cluster with V_1. Default value is 0.
	public double E_V = 0.0;

	// The binding energy of the cluster with I_1. Default value is 0.
	public double E_I = 0.0;

	// The migration energy with which the cluster diffuses. Default value is 0.
	public double E_m = 0.0;

	// The diffusion factor that sets the scale of the diffusion. Default value
	// is 0.
	public double D_0 = 0.0;

	/**
	 * This operation returns the cluster as a string of the form
	 * 
	 * nHe nV nI E_He E_V E_I E_migration D_0
	 */
	public String toString() {

		// Convert infinite energies
		String EHeAsString = (!Double.isInfinite(E_He)) ? String.valueOf(E_He)
				: "Infinity";
		String EVAsString = (!Double.isInfinite(E_V)) ? String.valueOf(E_V)
				: "Infinity";
		String EIAsString = (!Double.isInfinite(E_I)) ? String.valueOf(E_I)
				: "Infinity";
		String EMAsString = (!Double.isInfinite(E_I)) ? String.valueOf(E_m)
				: "Infinity";

		return nHe + " " + nV + " " + nI + " " + EHeAsString + " " + EVAsString
				+ " " + EIAsString + " " + EMAsString + " " + D_0;
	}

}
