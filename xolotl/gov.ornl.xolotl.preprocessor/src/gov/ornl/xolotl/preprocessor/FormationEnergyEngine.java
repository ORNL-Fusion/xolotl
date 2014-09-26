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
 * This class computes the formation energies for a cluster with a given
 * composition. All formation energies are in electron volts (eV).
 * 
 * This class looks in the user directory during construction and will override
 * the fit coefficients if it finds a csv file with one set of coefficients per
 * row and named fit.csv. There should be three rows with three coefficients per
 * row.
 * 
 * @author Jay Jay Billings
 * 
 */
public class FormationEnergyEngine {

	/**
	 * The set of Helium formation energies up to He_8 indexed by size. That is
	 * E_(f,He_1) = heFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] heFormationEnergies = { Double.POSITIVE_INFINITY,
			6.15, 11.44, 16.35, 21.0, 26.1, 30.24, 34.93, 38.80 };

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
	public FormationEnergyEngine() {

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
				System.err.println("FormationEnergyEngine Message:"
						+ " There is something wrong with your fit.csv file! "
						+ " It does not exist. Aborting.");
				e.printStackTrace();
			} catch (IOException e) {
				// Complain
				System.err.println("FormationEnergyEngine Message:"
						+ " There is something wrong with your fit.csv file! "
						+ " It cannot be read. Aborting.");
				e.printStackTrace();
			}
		}

		return;
	};

	/**
	 * This operation computes and returns the formation energy of a helium
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the Helium cluster.
	 * @return The formation energy.
	 */
	public double getHeFormationEnergy(int size) {

		double energy = Double.POSITIVE_INFINITY;

		if (size < 9 && size > 0)
			energy = heFormationEnergies[size];

		return energy;
	}

	/**
	 * This operation computes and returns the formation energy of a vacancy
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the vacancy cluster.
	 * @return The formation energy.
	 */
	public double getVFormationEnergy(int size) {

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
	 * This operation returns the formation energy for a cluster composed of only
	 * interstitial defects.
	 * 
	 * The values are taken from the iFormationEnergies array for a size smaller 
	 * or equal to 5, and computed with the formula
	 * 		E_f = 48 + 6 * (size - 6) eV
	 * starting at size = 6.
	 * 
	 * @param size
	 *            the number of interstitials in the cluster
	 * @return the formation energy
	 */
	public double getIFormationEnergy(int size) {
		// Initialize the energy
		double energy = Double.POSITIVE_INFINITY;

		// The formation energies from the table are valid for 0 < size < 6.
		if (size > 0 && size < 6) {
			energy = iFormationEnergies[size];
		} else if (size > 5)
			energy = 48.0 + 6.0 * (size - 6.0);

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
	public double getHeVFormationEnergy(int heSize, int vSize) {

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
			// Unscaling the energy
			energy = energy * 281.254;

		} else if ((vSize == 1 && heSize < heV1FormationEnergies.length)
				|| (vSize == 2 && heSize < heV2FormationEnergies.length)) {
			// Get the exact energy
			energy = (vSize == 1) ? heV1FormationEnergies[heSize]
					: heV2FormationEnergies[heSize];
		} 

		return energy;
	}

}
