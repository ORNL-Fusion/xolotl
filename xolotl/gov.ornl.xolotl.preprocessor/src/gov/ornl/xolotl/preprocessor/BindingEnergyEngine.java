/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * This class computes the binding energies for a cluster with a given
 * composition. All binding energies are in electron volts (eV).
 * 
 * This class looks in the user directory during construction and will override
 * the fit coefficients if it finds a csv file with one set of coefficients per
 * row and named fit.csv. There should be six rows with four coefficients per
 * row.
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
	private static double[] heFormationEnergies = { Double.POSITIVE_INFINITY,
			6.15, 11.44, 16.35, 21.0, 26.1, 30.24, 34.93, 38.80, 42.90 };

	/**
	 * The set of vacancy formation energies up to V_2 indexed by size. That is
	 * E_(f,V_1) = vFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] vFormationEnergies = { Double.POSITIVE_INFINITY,
			3.6, 7.25 };

	/**
	 * The set of interstitial formation energies up to I_6 indexed by size.
	 * That is E_(f,V_1) = iFormationEnergies[1]. The value at index zero is
	 * just padding to make the indexing easy.
	 */
	private static double[] iFormationEnergies = { Double.POSITIVE_INFINITY,
			10.0, 18.5, 27.0, 35.0, 42.5, 48.0 };

	// Coefficients for the Legendre polynomial fit
	// Coefficients for c_0 in the 2D E_f,HeV fit
	 private double[] c0Coefficients = { 8.11632e-01, 6.99475e-01,
	 -4.15369e-02 };
	 // Coefficients for c_1 in the 2D E_f,HeV fit
	 private double[] c1Coefficients = { 1.03498e+00, 9.53649e-01,
	 -2.91401e-02 };
	 // Coefficients for c_2 in the 2D E_f,HeV fit
	 private double[] c2Coefficients = { 3.63264e-01, 3.67345e-01,
	 -1.57605e-04 };

	/**
	 * The formation energies for He_xV_1. The entry at i = 0 is for a single
	 * vacancy (He_0V_1) and is there as a buffer. Like the formation energies,
	 * i = heSize.
	 */
	private double[] heV1FormationEnergies = { vFormationEnergies[1], 5.14166,
			8.20919, 11.5304, 14.8829, 18.6971, 22.2847, 26.3631, 30.1049,
			34.0081, 38.2069, 42.4217, 46.7378, 51.1551, 55.6738 };

	/**
	 * The formation energies for He_xV_2. The entry at i = 0 is for a
	 * di-vacancy (He_0V_2) and is there as a buffer. Like the formation
	 * energies, i = heSize.
	 */
	private double[] heV2FormationEnergies = { vFormationEnergies[2], 7.10098,
			8.39913, 9.41133, 11.8748, 14.8296, 17.7259, 20.7747, 23.7993,
			26.7984, 30.0626, 33.0385, 36.5173, 39.9406, 43.48, 46.8537,
			50.4484, 54.0879, 57.7939 };

	/**
	 * The constructor. If the file fit.csv exists in user.dir, this operation
	 * will try to load the coefficients from it. It closes the file when it is
	 * finished.
	 */
	public BindingEnergyEngine() {

		// Get the fit file
		File fitFile = new File("fit.csv");

		// Load the file if it exists
		if (fitFile.exists()) {
			try {
				// Create the readers
				FileReader fitFileReader = new FileReader(fitFile);
				BufferedReader fitReader = new BufferedReader(fitFileReader);

				// Read from the Legendre fit file
				String[] c0Strings = fitReader.readLine().split(",");
				String[] c1Strings = fitReader.readLine().split(",");
				String[] c2Strings = fitReader.readLine().split(",");

				// Convert the c0 coefficients
				c0Coefficients[0] = Double.valueOf(c0Strings[0]);
				c0Coefficients[1] = Double.valueOf(c0Strings[1]);
				c0Coefficients[2] = Double.valueOf(c0Strings[2]);
				// Convert the c1 coefficients
				c1Coefficients[0] = Double.valueOf(c1Strings[0]);
				c1Coefficients[1] = Double.valueOf(c1Strings[1]);
				c1Coefficients[2] = Double.valueOf(c1Strings[2]);
				// Convert the c2 coefficients
				c2Coefficients[0] = Double.valueOf(c2Strings[0]);
				c2Coefficients[1] = Double.valueOf(c2Strings[1]);
				c2Coefficients[2] = Double.valueOf(c2Strings[2]);

				// Close the reader
				fitReader.close();
			} catch (FileNotFoundException e) {
				// Complain
				System.err.println("BindingEnergyEngine Message:"
						+ " There is something wrong with your fit.csv file! "
						+ " It does not exist. Aborting.");
				e.printStackTrace();
			} catch (IOException e) {
				// Complain
				System.err.println("BindingEnergyEngine Message:"
						+ " There is something wrong with your fit.csv file! "
						+ " It cannot be read. Aborting.");
				e.printStackTrace();
			}
		}

		return;
	};

	/**
	 * This operation computes and returns the formation energy of a Helium
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the Helium cluster.
	 * @return The formation energy or the proper type of infinite to satisfy
	 *         the binding energy equation.
	 */
	private double getHeFormationEnergy(int size) {

		double energy = Double.POSITIVE_INFINITY;

		if (size <= 9 && size > 0)
			energy = heFormationEnergies[size];
		else if (size > 9)
			energy = Double.NEGATIVE_INFINITY;

		return energy;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of only
	 * Helium that dissociates as
	 * 
	 * He_x --> He_(x-1) + He_1
	 * 
	 * For size > 8 or size < 2 it returns Double.POSITIVE_INFINITY.
	 * 
	 * @param size
	 *            the number of Helium atoms in the dissociating cluster
	 * @return the binding energy
	 */
	public double getHeBindingEnergy(int size) {
		return getHeFormationEnergy(size - 1) + heFormationEnergies[1]
				- getHeFormationEnergy(size);
	}

	/**
	 * This operation computes and returns the formation energy of a vacancy
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the vacancy cluster.
	 * @return The formation energy or the proper type of infinite to satisfy
	 *         the binding energy equation.
	 */
	private double getVFormationEnergy(int size) {

		double energy = Double.POSITIVE_INFINITY;

		if (size < 3 && size > 0)
			energy = vFormationEnergies[size];
		else if (size >= 3) {
			// Get the energy for He_xV_y
			energy = getHeVFormationEnergy(0, size);
		}

		return energy;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of only
	 * atomic vacancies that dissociates as
	 * 
	 * V_x --> V_(x-1) + V_1
	 * 
	 * It is valid for all sizes of V > 0. For size > 2 it uses the same fit as
	 * that used for HeV clusters.
	 * 
	 * @param size
	 *            the number of atomic vacancies in the dissociating cluster
	 * @return the binding energy
	 */
	public double getVBindingEnergy(int size) {
		return getVFormationEnergy(size - 1) + vFormationEnergies[1]
				- getVFormationEnergy(size);
	}

	/**
	 * This operation returns the binding energy for a cluster composed of only
	 * interstitial defects that dissociates as
	 * 
	 * I_x --> I_(x-1) + I_1
	 * 
	 * For size < 2 it returns Double.POSITIVE_INFINITY. It is valid for all
	 * other sizes and computes the formation energy using
	 * 
	 * E_f(size) = 48 + 6*(size - 6)
	 * 
	 * @param size
	 *            the number of interstitials in the cluster
	 * @return the binding energy
	 */
	public double getIBindingEnergy(int size) {

		double energy = Double.POSITIVE_INFINITY;

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
	 * This operation computes the Legendre polynomial, P_n (x), of degree n.
	 * 
	 * The Legendre polynomials of degree 0 and 1 are P_0(x)=1.0 and P_1(x)=x,
	 * respectively. With these conditions, the Legendre polynomials satisfy the
	 * following recurrence relation: 
	 * (n+1)*P_(n+1)(x) = (2n+1)*x*P_n(x) - n*P_(n-1)(x)
	 * 
	 * @param x
	 *            The x value of the function
	 * @param order
	 *            The order of the polynomial
	 */
	private double legendrePolynomial(double x, int degree) {

		if (degree == 0)
			return 1.0;
		else if (degree == 1)
			return x;

		double Pn2 = 1.0, Pn1 = x, Pn = 0.0;

		for (int n = 1; n < degree; n++) {
			Pn = (((2.0 * (double) n + 1.0) * x * Pn1) - ((double) n * Pn2))
					/ ((double) n + 1.0);
			Pn2 = Pn1;
			Pn1 = Pn;
		}

		return Pn;
	}

	/**
	 * This operation computes the 2nd order Legendre polynomial expansion
	 * 
	 * f(x) = c0*P_0(x) + c1*P_1(x) + c2*P_2(x) = c0 + c1*x + c2*P_2(x)
	 * 
	 * for a coefficient set {c0,c1,c2}.
	 * 
	 * @param x
	 *            The x value of the function.
	 * @param coeffs
	 *            The coefficients array.
	 */
	private double compute2ndOrderLegendrePCE(double x, double[] coeffs) {

		double value = 0.0;

		// Compute the value
		value = coeffs[0] + coeffs[1] * x + coeffs[2]
				* legendrePolynomial(x, 2);

		return value;
	}

	/**
	 * This operation returns the formation energy for an HeV cluster.
	 * 
	 * The formation energies are determined using the following 2D,
	 * second order, Legendre polynomial fit
	 * 
	 * f(y,x/y) = [c_0 + c_1*y + c_2*P_2(y)] + [c_3 + c_4*y + c_5*P_2(y)]*(x/y)
	 * 			 + [c_6 + c_7*y + c_8*P_2(y)]*P_2(x/y)
	 * 
	 * where y = vSize and x/y = heSize/vSize
	 * 
	 * @param heSize
	 *            The number of Helium atoms in the cluster
	 * @param vSize
	 *            The number of vacancies in the cluster
	 * @return The formation energy
	 */
	private double getHeVFormationEnergy(int heSize, int vSize) {

		double energy = Double.NEGATIVE_INFINITY;
		// The following coefficients are computed using the above and are used
		// to evaluate the full function f(x,y).
		double[] coefficients = { 0.0, 0.0, 0.0 };

		// Check to see if the vacancy size is large enough that the energy can
		// be computed from the fit or if it is so small that the exact values
		// must be used instead.
		if (vSize > 2) {
			// Get the Vacancy number
			double y = 2.0 * (((double) vSize - 1.0) / 43.0) - 1.0;
			// Get the coefficients
			coefficients[0] = compute2ndOrderLegendrePCE(y,
					c0Coefficients);
			coefficients[1] = compute2ndOrderLegendrePCE(y,
					c1Coefficients);
			coefficients[2] = compute2ndOrderLegendrePCE(y,
					c2Coefficients);
			// Get the He/V ratio
			double x = 2.0 * (((double) heSize / (double) vSize) / 9.0) - 1.0;
			// Get the energy for He_xV_y
			energy = compute2ndOrderLegendrePCE(x, coefficients);

			energy = energy * 281.254;

		} else if ((vSize == 1 && heSize < heV1FormationEnergies.length)
				|| (vSize == 2 && heSize < heV2FormationEnergies.length)) {
			// Get the exact energy
			energy = (vSize == 1) ? heV1FormationEnergies[heSize]
					: heV2FormationEnergies[heSize];
		}

		return energy;
	}

	/**
	 * This operation returns the binding energy for a cluster composed of
	 * Helium and vacancies that dissociates as
	 * 
	 * (He_x)(V_y) --> (He_x-1)(V_y) + He_1
	 * 
	 * It uses a 2D Legendre polynomial fit of the formation energies as a
	 * function of vsize and heSize/vSize to compute the formation energy 
	 * and then the binding energy.  The fit is of the form
	 * 
	 * f(y,x/y) = c_0 + c_1*(x/y) + c_2*P_2(x/y)
	 * 
	 * where each coefficient c_i is determined as follows,
	 * 
	 * c_i = d_0 + d_1*y + d_2*P_2(y)
	 * 
	 * Hence,
	 * 
	 * f(y,x/y) = [c_0 + c_1*y + c_2*P_2(y)] + [c_3 + c_4*y + c_5*P_2(y)]*(x/y)
	 * 			 + [c_6 + c_7*y + c_8*P_2(y)]*P_2(x/y)
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
		double heVE_f = 0.0, heMinusOneVE_f = 0.0, energy = Double.POSITIVE_INFINITY;

		// Check the bounds and delegate if necessary. Start with the normal
		// case where both heSize and vSize are equal to or greater than 1.
		if (heSize > 0 && vSize > 0) {
			// Get the energy for He_xV_y
			heVE_f = getHeVFormationEnergy(heSize, vSize);
			// Most of the time the clusters will have more than one Helium
			// atom.
			if (heSize > 1) {
				// Get the energy for He_(x-1)V_y
				heMinusOneVE_f = getHeVFormationEnergy(heSize - 1, vSize);
				// Compute the energy
				energy = heMinusOneVE_f + heFormationEnergies[1] - heVE_f;
			} else {
				// heSize == 1 is a special case where the cluster breaks apart
				// into the vacancy clusters plus single helium.
				energy = heFormationEnergies[1] + getVFormationEnergy(vSize)
						- heVE_f;
			}
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
	 * It uses a 2D Legendre polynomial fit of the formation energies as a
	 * function of vsize and heSize/vSize to compute the formation energy 
	 * and then the binding energy.  The fit is of the form
	 * 
	 * f(y,x/y) = c_0 + c_1*(x/y) + c_2*P_2(x/y)
	 * 
	 * where each coefficient c_i is determined as follows,
	 * 
	 * c_i = d_0 + d_1*y + d_2*P_2(y)
	 * 
	 * Hence,
	 * 
	 * f(y,x/y) = [c_0 + c_1*y + c_2*P_2(y)] + [c_3 + c_4*y + c_5*P_2(y)]*(x/y)
	 * 			 + [c_6 + c_7*y + c_8*P_2(y)]*P_2(x/y)
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
		double heVE_f = 0.0, heVMinusOneE_f = 0.0, energy = Double.POSITIVE_INFINITY;

		// Check the bounds and delegate if necessary
		if (heSize > 0 && vSize > 0) {
			// Get the energy for He_xV_y
			heVE_f = getHeVFormationEnergy(heSize, vSize);
			// Most of the time the clusters will have more than one Helium
			// atom.
			if (vSize > 1) {
				// Get the energy for He_(x-1)V_y. This accounts for the case
				// where V == 1 by getting the He binding energy instead of
				// querying
				// the fit, which results in a NaN in that case. This is a
				// pretty
				// special case in general and needs to be investigated further.
				heVMinusOneE_f = (vSize == 1) ? getHeFormationEnergy(heSize)
						: getHeVFormationEnergy(heSize, vSize - 1);
				// Compute the energy
				energy = heVMinusOneE_f + vFormationEnergies[1] - heVE_f;
			} else {
				// vSize == 1 is a special case where the cluster breaks apart
				// into the vacancy clusters plus single helium.
				energy = getHeFormationEnergy(heSize) + vFormationEnergies[1]
						- heVE_f;
			}
		} else if (heSize == 0) {
			// Deal with vacancies
			energy = getVBindingEnergy(vSize);
		}

		return energy;
	}

	// Add binding energy routines for HeI below.

}
