/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

/**
 * This class computes the binding energies for a cluster with a given
 * composition. All binding energies are in electron volts (eV).
 * 
 * @author Jay Jay Billings
 * 
 */
public class BindingEnergyEngine {

	/**
	 * The set of Helium formation energies up to He_8 indexed by size. That is
	 * E_(f,He_1) = heFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] heFormationEnergies = { Double.MAX_VALUE, 6.15,
			11.44, 16.35, 21.0, 26.1, 30.24, 34.93, 38.80 };

	/**
	 * The set of vacancy formation energies up to V_2 indexed by size. That is
	 * E_(f,V_1) = vFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] vFormationEnergies = { Double.MAX_VALUE, 3.60, 7.25 };

	/**
	 * The set of interstitial formation energies up to I_6 indexed by size.
	 * That is E_(f,V_1) = iFormationEnergies[1]. The value at index zero is
	 * just padding to make the indexing easy.
	 */
	private static double[] iFormationEnergies = { Double.MAX_VALUE, 10.0,
			18.5, 27.0, 35.0, 42.5, 48.0 };

	// Coefficients for c_0 in the 2D E_f,HeV fit
	private double[] c0Coefficients = { 2.0918047239, 2.424457391, -0.0330724815,
			0.0002688409 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	private double[] c1Coefficients = { 1.3409500982, -0.8075623963, 0.0378738095,
			-0.0004038271 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	private double[] c2Coefficients = { 1.0115538203, 0.9313017002, -0.0612676135,
			0.0006779278 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	private double[] c3Coefficients = { -0.3332680916, -0.2007444968, 0.0320255177,
			-0.0003860637 };
	// Coefficients for c_4 in the 2D E_f,HeV fit
	private double[] c4Coefficients = { 0.0145638563, 0.0466697556, -0.007875649,
			0.0001003462 };
	// Coefficients for c_5 in the 2D E_f,HeV fit
	private double[] c5Coefficients = { 0.0034776952, -0.0050987442, 0.0007179606,
			-9.40486423253e-6 };
	
	/**
	 * This operation returns the binding energy for a cluster composed of only
	 * Helium that dissociates as
	 * 
	 * He_x --> He_(x-1) + He_1
	 * 
	 * For size > 8 or size < 2 it returns Double.MAX_VALUE.
	 * 
	 * @param size
	 *            the number of Helium atoms in the dissociating cluster
	 * @return the binding energy
	 */
	public double getHeBindingEnergy(int size) {

		double energy = (size < 9 && size > 1) ? heFormationEnergies[size - 1]
				+ heFormationEnergies[1] - heFormationEnergies[size]
				: Double.MAX_VALUE;

		return energy;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of only
	 * atomic vacancies that dissociates as
	 * 
	 * V_x --> V_(x-1) + V_1
	 * 
	 * For size != 2 it returns Double.MAX_VALUE.
	 * 
	 * @param size
	 *            the number of atomic vacancies in the dissociating cluster
	 * @return the binding energy
	 */
	public double getVBindingEnergy(int size) {

		double energy = (size == 2) ? vFormationEnergies[size - 1]
				+ vFormationEnergies[1] - vFormationEnergies[size]
				: Double.MAX_VALUE;

		return energy;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of only
	 * interstitial defects that dissociates as
	 * 
	 * I_x --> I_(x-1) + I_1
	 * 
	 * For size < 2 it returns Double.MAX_VALUE. It is valid for all other sizes
	 * and computes the formation energy using
	 * 
	 * E_f(size) = 48 + 6*(size - 6)
	 * 
	 * @param size
	 *            the number of interstitials in the cluster
	 * @return the binding energy
	 */
	public double getIBindingEnergy(int size) {

		double energy = Double.MAX_VALUE;

		// The formation energies from the table are valid for 2 <= size < 6.
		if (size >= 2 && size < 6) {
			energy = iFormationEnergies[size - 1] + iFormationEnergies[1]
					- iFormationEnergies[size];
		} else if (size == 6) {
			// The formation energy for I_6 is 48 based on the formula.
			energy = 48.0 + iFormationEnergies[1] - iFormationEnergies[5];
		} else if (size > 6)
			// All larger clusters have the same binding energy based on
			// the formula.
			energy = 16.0;

		return energy;
	}

	/**
	 * This operation computes the function
	 * 
	 * f(x) = a + b*x + c*x^2 + d*x^3
	 * 
	 * for a coefficient set {a,b,c,d}.
	 * 
	 * This function is called to compute coefficients c_i for the vacancy size
	 * dimension of the HeV binding energy fit.
	 * 
	 * @param x
	 *            The x value of the function.
	 * @param coeffs
	 *            The coefficients array.
	 */
	private double compute3rdOrderPolynomial(double x, double[] coeffs) {

		// Local Declarations
		double value = 0.0;

		// Compute the value
		value = coeffs[0] + coeffs[1] * x + coeffs[2] * x * x + coeffs[3] * x
				* x * x;

		return value;
	}

	/**
	 * This operation computes the function
	 * 
	 * f(x) = a + b*x + c*x^2 + d*x^3 + e*x^4 + f*x^5
	 * 
	 * for a coefficient set {a,b,c,d,e,f}.
	 * 
	 * This function is called to compute the final formation energy for the HeV
	 * formation energy fit f(x,y).
	 * 
	 * @param x
	 *            The x value of the function.
	 * @param coeffs
	 *            The coefficients array.
	 */
	private double compute5thOrderPolynomial(double x, double[] coeffs) {

		// Local Declarations
		double value = 0.0;

		// Compute the value
		value = coeffs[0] + coeffs[1] * x + coeffs[2] * x * x + coeffs[3] * x
				* x * x + coeffs[4] * x * x * x * x + coeffs[5] * x * x * x * x
				* x;

		return value;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of
	 * Helium and vacancies that dissociates as
	 * 
	 * (He_x)(V_y) --> (He_x-1)(V_y) + He_1
	 * 
	 * It uses a 2D fit of the formation energies as a function of vSize and
	 * heSize/vSize to compute the formation energy and then the binding energy.
	 * The fit is of the form
	 * 
	 * f(y,y/x) = c_0(y) + c_1(y)*(y/x) + c_2(y)*(y/x)^2 + c_3(y)*(y/x)^3 +
	 * c_4(y)*(y/x)^4 + c_5(y)*(y/x)^5 + c_6(y)*(y/x)^6
	 * 
	 * where each function c_i(y) is
	 * 
	 * c_i(y) = d_0 + d_1*y + d_2*y^2 + d_3*y_3 + d_4*y_4
	 * 
	 * If heSize or vSize is equal to zero, this function will return the value
	 * of the single species routines for He and V respectively.
	 * 
	 * @param heSize
	 *            the number of Helium atoms in the cluster
	 * @param vSize
	 *            the number of vacancies in the cluster
	 * @return the binding energy
	 */
	public double getHeVtoHeBindingEnergy(int heSize, int vSize) {

		// Local Declarations
		double heVE_f = 0.0, heMinusOneVE_f = 0.0, energy = 0.0;
		// The following coefficients are computed using the above and are used
		// to evaluate the full function f(x,y).
		double[] coefficients = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

		// Check the bounds and delegate if necessary
		if (heSize > 0 && vSize > 0) {
			// Get the coefficients
			coefficients[0] = compute3rdOrderPolynomial((double) vSize, c0Coefficients);
			coefficients[1] = compute3rdOrderPolynomial((double) vSize, c1Coefficients);
			coefficients[2] = compute3rdOrderPolynomial((double) vSize, c2Coefficients);
			coefficients[3] = compute3rdOrderPolynomial((double) vSize, c3Coefficients);
			coefficients[4] = compute3rdOrderPolynomial((double) vSize, c4Coefficients);
			coefficients[5] = compute3rdOrderPolynomial((double) vSize, c5Coefficients);
			// Get the energy for He_xV_y
			heVE_f = compute5thOrderPolynomial((double) heSize / (double) vSize, coefficients);
			// Get the energy for He_(x-1)V_y
			heMinusOneVE_f = compute5thOrderPolynomial((double) (heSize - 1) / (double) vSize,
					coefficients);
			// Compute the energy
			energy = heMinusOneVE_f + heFormationEnergies[1] - heVE_f;
		} else if (heSize == 0) {
			// Deal with vacancies
			energy = getVBindingEnergy(vSize);
		} else if (vSize == 0) {
			// Deal with Helium
			energy = getHeBindingEnergy(heSize);
		}

		return energy;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of
	 * Helium and vacancies that dissociates as
	 * 
	 * (He_x)(V_y) --> (He_x)(V_y-1) + V_1
	 * 
	 * It uses a 2D fit of the formation energies as a function of vSize and
	 * heSize/vSize to compute the formation energy and then the binding energy.
	 * 
	 * f(y,y/x) = c_0(y) + c_1(y)*(y/x) + c_2(y)*(y/x)^2 + c_3(y)*(y/x)^3 +
	 * c_4(y)*(y/x)^4 + c_5(y)*(y/x)^5 + c_6(y)*(y/x)^6
	 * 
	 * where each function c_i(y) is
	 * 
	 * c_i(y) = d_0 + d_1*y + d_2*y^2 + d_3*y_3 + d_4*y_4
	 * 
	 * If heSize or vSize is equal to zero, this function will return the value
	 * of the single species routines for He and V respectively.
	 * 
	 * @param heSize
	 *            the number of Helium atoms in the cluster
	 * @param vSize
	 *            the number of vacancies in the cluster
	 * @return the binding energy
	 */
	public double getHeVtoVBindingEnergy(int heSize, int vSize) {
		// Local Declarations
		double heVE_f = 0.0, heVMinusOneE_f = 0.0, energy = 0.0;
		// The following coefficients are computed using the above and are used
		// to evaluate the full function f(x,y).
		double[] coefficients = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

		// Check the bounds and delegate if necessary
		if (heSize > 0 && vSize > 0) {
			// Get the coefficients for He_xV_y
			coefficients[0] = compute3rdOrderPolynomial((double) vSize, c0Coefficients);
			coefficients[1] = compute3rdOrderPolynomial((double) vSize, c1Coefficients);
			coefficients[2] = compute3rdOrderPolynomial((double) vSize, c2Coefficients);
			coefficients[3] = compute3rdOrderPolynomial((double) vSize, c3Coefficients);
			coefficients[4] = compute3rdOrderPolynomial((double) vSize, c4Coefficients);
			coefficients[5] = compute3rdOrderPolynomial((double) vSize, c5Coefficients);
			// Get the energy for He_xV_y
			heVE_f = compute5thOrderPolynomial((double) heSize / (double) vSize, coefficients);
			// Recompute the coefficients to compute He_xV_(y-1)
			double vMinus1 = (double) (vSize-1);
			coefficients[0] = compute3rdOrderPolynomial(vMinus1, c0Coefficients);
			coefficients[1] = compute3rdOrderPolynomial(vMinus1, c1Coefficients);
			coefficients[2] = compute3rdOrderPolynomial(vMinus1, c2Coefficients);
			coefficients[3] = compute3rdOrderPolynomial(vMinus1, c3Coefficients);
			coefficients[4] = compute3rdOrderPolynomial(vMinus1, c4Coefficients);
			coefficients[5] = compute3rdOrderPolynomial(vMinus1, c5Coefficients);
			// Get the energy for He_(x-1)V_y
			heVMinusOneE_f = compute5thOrderPolynomial((double) heSize / vMinus1,
					coefficients);
			// Compute the energy
			energy = heVMinusOneE_f + vFormationEnergies[1] - heVE_f;
		} else if (heSize == 0) {
			// Deal with vacancies
			energy = getVBindingEnergy(vSize);
		} else if (vSize == 0) {
			// Deal with Helium
			energy = getHeBindingEnergy(heSize);
		}

		return energy;
	}

	// Add binding energy routines for HeI below.

}
